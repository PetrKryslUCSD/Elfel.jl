module heat_poisson_t3

using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
# using Profile
using MeshCore: retrieve, nrelations, nentities
using MeshSteward: T3block
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: connectedv, geometry
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T3, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, numberdofs!, setebc!, nunknowns, doftype
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: asstolma!, lma, asstolva!, lva
using Elfel.QPIterators: QPIterator, bfun, bfungradpar, weight, jacjac
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler

A = 1.0 # length of the side of the square
kappa =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = T3block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function assembleK(fesp, kappa)
    function integrateK!(ass, geom, elit, qpit, kappa)
        nedof = ndofsperel(elit)
        for el in elit
            for qp in qpit
                gradNparams = bfungradpar(qp)
                Jac, J = jacjac(qp, geom, elnodes(el))
                JxW = J * weight(qp)
                invJac = inv(Jac)
                for j in 1:nedof
                    gradNj = gradNparams[j] * invJac
                    for i in 1:nedof
                        gradNi = gradNparams[i] * invJac
                        v = dot(gradNi, gradNj) * kappa * JxW 
                        asstolma!(el, i, j, v)
                    end
                end
            end
            assemble!(ass, lma(el)...)
        end
        return ass
    end

    elit = FEIterator(fesp)
    qpit = QPIterator(fesp.fe, (kind = :default,))
    geom = geometry(fesp.mesh)
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, ndofs(fesp), ndofs(fesp))
    @time integrateK!(ass, geom, elit, qpit, kappa)
    return finish!(ass)
end

function assembleF(fesp, Q)
    function integrateF!(ass, geom, elit, qpit, kappa)
        nedof = ndofsperel(elit)
        for el in elit
            for qp in qpit
                Jac, J = jacjac(qp, geom, elnodes(el))
                JxW = J * weight(qp)
                N = bfun(qp)
                for i in 1:nedof
                    asstolva!(el, i, N[i] * Q * JxW)
                end
            end
            assemble!(ass, lva(el)...)
        end
        return ass
    end

    elit = FEIterator(fesp)
    qpit = QPIterator(fesp.fe, (kind = :default,))
    geom = geometry(fesp.mesh)
    ass = SysvecAssembler(0.0)
    start!(ass, ndofs(fesp))
    @time integrateF!(ass, geom, elit, qpit, Q)
    return finish!(ass)
end

function solve!(T, K, F, nu)
    @time KT = K * T
    @time T[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function run()
    mesh = genmesh()
    fesp = FESpace(Float64, FEH1_T3(1), mesh)
    bir = boundary(mesh);
    vl = connectedv(bir);
    locs = geometry(mesh)
    for i in vl
        setebc!(fesp, 0, i, 1, tempf(locs[i]...))
    end
    numberdofs!(fesp)
    @show nunknowns(fesp)
    K = assembleK(fesp, kappa)
    F = assembleF(fesp, Q)
    T = fill(0.0, ndofs(fesp))
    gathersysvec!(T, fesp)
    solve!(T, K, F, nunknowns(fesp))
    scattersysvec!(fesp, T)
    makeattribute(fesp, "T", 1)
    vtkwrite("heat_poisson_t3-T", baseincrel(mesh), ["T"])
end

end

heat_poisson_t3.run()
# heat_poisson_t3.run()
