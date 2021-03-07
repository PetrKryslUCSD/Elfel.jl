"""
    q4

Compute the solution of the Poisson equation of heat conduction with a nonzero
heat source. Quadrilateral four-node elements are used.
"""
module q4

using LinearAlgebra
using SparseArrays
using MeshCore.Exports
using MeshCore: @_check
using MeshSteward.Exports
using Elfel.Exports
using UnicodePlots

A = 1.0 # length of the side of the square
kappa =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = Q4block(A, A, N, N)
    mesh = Mesh()
    attach!(mesh, conn)
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

    @time integrate!(am, av, geom, elit, qpit, kappa, Q)

    return finish!(am), finish!(av)
end

function solve!(T, K, F, nu)
    @time KT = K * T
    @time T[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function checkcorrectness(fesp)
    geom = geometry(fesp.mesh)
    ir = baseincrel(fesp.mesh)
    T = attribute(ir.right, "T")
    std = 0.0
    for i in 1:length(T)
        std += abs(T[i][1] - tempf(geom[i]...))
    end
    @_check (std / length(T)) <= 1.0e-9
end

function run()
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
    @show nunknowns(fesp)
    K, F = assembleKF(fesp, kappa, Q)
    T = fill(0.0, ndofs(fesp))
    gathersysvec!(T, fesp)
    solve!(T, K, F, nunknowns(fesp))
    scattersysvec!(fesp, T)
    makeattribute(fesp, "T", 1)
    checkcorrectness(fesp)
    vtkwrite("q4-T", baseincrel(mesh), [(name = "T",)])
end

end

@time q4.run()
# q4.run()
