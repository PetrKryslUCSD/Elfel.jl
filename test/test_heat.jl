module mt_heat_poisson_t3

using Test
using LinearAlgebra
using MeshCore: retrieve, nrelations, nentities, @_check, attribute
using MeshSteward: T3block
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: connectedv, geometry
using MeshSteward: vtkwrite
using Elfel.RefShapes: RefShapeTriangle, manifdim, manifdimv
using Elfel.FElements: FEH1_T3, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: jacjac
using Elfel.QPIterators: QPIterator, bfun, bfungradpar, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!

A = 1.0 # length of the side of the square
kappa =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
N = 4;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = T3block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function assembleKF(fesp, kappa, Q)
    function integrate!(am, av, geom, elit, qpit, kappa, Q)
        nedof = ndofsperel(elit)
        iterate(qpit, 1) # Why do I need to do this?
        ke = LocalMatrixAssembler(nedof, nedof, 0.0)
        fe = LocalVectorAssembler(nedof, 0.0)
        for el in elit
            init!(ke, eldofs(el), eldofs(el))
            init!(fe, eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        ke[i, j] += dot(gradN[i], gradN[j]) * (kappa * JxW)
                    end
                    fe[j] += N[j] * Q * JxW
                end
            end
            assemble!(am, ke)
            assemble!(av, fe)
        end
        return am, av
    end

    elit = FEIterator(fesp)
    qpit = QPIterator(fesp, (kind = :default,))
    geom = geometry(fesp.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(fesp), ndofs(fesp))
    av = start!(SysvecAssembler(0.0), ndofs(fesp))

    integrate!(am, av, geom, elit, qpit, kappa, Q)

    return finish!(am), finish!(av)
end

function solve!(T, K, F, nu)
    KT = K * T
    T[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function checkcorrectness(fesp)
    geom = geometry(fesp.mesh)
    ir = baseincrel(fesp.mesh)
    T = attribute(ir.right, "T")
    std = 0.0
    for i in 1:length(T)
        std += abs(T[i][1] - tempf(geom[i]...))
    end
    @test (std / length(T)) <= 1.0e-9
end

function test()
    mesh = genmesh()
    fesp = FESpace(Float64, mesh, FEH1_T3())
    bir = boundary(mesh);
    vl = connectedv(bir);
    locs = geometry(mesh)
    for i in vl
        setebc!(fesp, 0, i, 1, tempf(locs[i]...))
    end
    numberfreedofs!(fesp)
    numberdatadofs!(fesp)
    # @show nunknowns(fesp)
    K, F = assembleKF(fesp, kappa, Q)
    T = fill(0.0, ndofs(fesp))
    gathersysvec!(T, fesp)
    solve!(T, K, F, nunknowns(fesp))
    scattersysvec!(fesp, T)
    makeattribute(fesp, "T", 1)
    checkcorrectness(fesp)
    vtkwrite("heat_poisson_t3-T", baseincrel(mesh), [(name = "T",)])
    try rm("heat_poisson_t3-T.vtu"); catch end
    @test isapprox(T, [1.1875, 1.3749999999999998, 1.6874999999999998, 1.5624999999999998, 1.7499999999999998, 2.0625, 2.1875, 2.375, 2.6875, 1.0, 1.0625, 1.25, 1.5625, 2.0, 1.125, 2.125, 1.5, 2.5, 2.125, 3.125, 3.0, 3.0625, 3.25, 3.5625, 4.0])
end

end
using .mt_heat_poisson_t3
mt_heat_poisson_t3.test()

module m_heat_poisson_q4
using Test
using LinearAlgebra
using MeshCore: retrieve, nrelations, nentities, attribute, @_check
using MeshSteward: Q4block
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: connectedv, geometry
using MeshSteward: vtkwrite
using Elfel.RefShapes: RefShapeTriangle, manifdim, manifdimv
using Elfel.FElements: FEH1_Q4, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: jacjac
using Elfel.QPIterators: QPIterator, bfun, bfungradpar, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!

A = 1.0 # length of the side of the square
kappa =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
N = 100;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = Q4block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function assembleKF(fesp, kappa, Q)
    function integrate!(am, av, geom, elit, qpit, kappa, Q)
        nedof = ndofsperel(elit)
        ke = LocalMatrixAssembler(nedof, nedof, 0.0)
        fe = LocalVectorAssembler(nedof, 0.0)
        for el in elit
            init!(ke, eldofs(el), eldofs(el))
            init!(fe, eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        ke[i, j] += dot(gradN[i], gradN[j]) * (kappa * JxW)
                    end
                    fe[j] += N[j] * Q * JxW
                end
            end
            assemble!(am, ke)
            assemble!(av, fe)
        end
        return am, av
    end

    elit = FEIterator(fesp)
    qpit = QPIterator(fesp, (kind = :Gauss, order = 2))
    geom = geometry(fesp.mesh)
    am = start!(SysmatAssemblerSparse(0.0), ndofs(fesp), ndofs(fesp))
    av = start!(SysvecAssembler(0.0), ndofs(fesp))

    integrate!(am, av, geom, elit, qpit, kappa, Q)

    return finish!(am), finish!(av)
end

function solve!(T, K, F, nu)
    KT = K * T
    T[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function checkcorrectness(fesp)
    geom = geometry(fesp.mesh)
    ir = baseincrel(fesp.mesh)
    T = attribute(ir.right, "T")
    std = 0.0
    for i in 1:length(T)
        std += abs(T[i][1] - tempf(geom[i]...))
    end
    @test (std / length(T)) <= 1.0e-9
end

function test()
    mesh = genmesh()
    fesp = FESpace(Float64, mesh, FEH1_Q4())
    bir = boundary(mesh);
    vl = connectedv(bir);
    locs = geometry(mesh)
    for i in vl
        setebc!(fesp, 0, i, 1, tempf(locs[i]...))
    end
    numberfreedofs!(fesp)
    numberdatadofs!(fesp)
    @test nunknowns(fesp) == 9801
    K, F = assembleKF(fesp, kappa, Q)
    T = fill(0.0, ndofs(fesp))
    gathersysvec!(T, fesp)
    solve!(T, K, F, nunknowns(fesp))
    scattersysvec!(fesp, T)
    makeattribute(fesp, "T", 1)
    checkcorrectness(fesp)
    vtkwrite("m_heat_poisson_q4-T", baseincrel(mesh), [(name = "T",)])
    try rm("m_heat_poisson_q4-T.vtu"); catch end
end

end
using .m_heat_poisson_q4
m_heat_poisson_q4.test()
