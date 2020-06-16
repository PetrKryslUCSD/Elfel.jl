module elast_stretch_t6

using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
using StaticArrays
# using Profile
using MeshCore: retrieve, nrelations, nentities
using MeshSteward: T6block
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: vselect, geometry
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T6, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, numberdofs!, setebc!, nunknowns, doftype
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!, edofcompnt
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: jacjac
using Elfel.QPIterators: QPIterator, bfun, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!

E = 1.0;
nu = 1.0/3;
D = SMatrix{3, 3}(E / (1 - nu^2) * [1 nu 0
                                    nu 1 0
                                    0 0 (1 - nu) / 2])
A = 1.0 # length of the side of the square
N = 10;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = T6block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function assembleK(fesp, D)
    function integrateK!(ass, geom, elit, qpit, D)
        B = (g, k) -> k == 1 ? SVector{3}((g[1], 0, g[2])) : SVector{3}((0, g[2], g[1]))
        c = edofcompnt(elit.fesp)
        nedof = ndofsperel(elit)
        ke = LocalMatrixAssembler(nedof, nedof, 0.0)
        for el in elit
            init!(ke, eldofs(el), eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                for j in 1:nedof
                    DBj = D * B(gradN[j], c[j])
                    for i in 1:nedof
                        Bi = B(gradN[i], c[i])
                        ke[i, j] += dot(DBj, Bi) * JxW
                    end
                end
            end
            assemble!(ass, ke)
        end
        return ass
    end

    elit = FEIterator(fesp)
    qpit = QPIterator(fesp, (kind = :default, npts = 3,))
    geom = geometry(fesp.mesh)
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, ndofs(fesp), ndofs(fesp))
    @time integrateK!(ass, geom, elit, qpit, D)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    @time KT = K * U
    @time U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function run()
    mesh = genmesh()
    fesp = FESpace(Float64, mesh, FEH1_T6(), 2)
    locs = geometry(mesh)
    inflate = A / N / 100
    box = [0.0 0.0 0.0 A]
    vl = vselect(locs; box = box, inflate = inflate)
    for i in vl
        setebc!(fesp, 0, i, 1, 0.0)
        setebc!(fesp, 0, i, 2, 0.0)
    end
    box = [A A 0.0 A]
    vl = vselect(locs; box = box, inflate = inflate)
    for i in vl
        setebc!(fesp, 0, i, 1, A / 10)
        setebc!(fesp, 0, i, 2, 0.0)
    end
    numberdofs!(fesp)
    @show nunknowns(fesp)
    K = assembleK(fesp, D)
    U = fill(0.0, ndofs(fesp))
    gathersysvec!(U, fesp)
    F = fill(0.0, ndofs(fesp))
    solve!(U, K, F, nunknowns(fesp))
    scattersysvec!(fesp, U)
    makeattribute(fesp, "U", 1:2)
    vtkwrite("elast_stretch_t6-U", baseincrel(mesh), [(name = "U", allxyz = true)])
end

end

elast_stretch_t6.run()
