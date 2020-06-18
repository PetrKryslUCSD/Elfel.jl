module stokes_driven_tht6

using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities
using MeshSteward: T6block, T6toT3
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: vselect, geometry, summary
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T6, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!, edofcompnt
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: jacjac
using Elfel.AggregateFEIterators: AggregateFEIterator, iterator
using Elfel.QPIterators: QPIterator, bfun, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!

mu = 0.25 # dynamic viscosity
A = 1.0 # length of the side of the square
N = 4;# number of subdivisions along the sides of the square domain

function genmesh()
    # Taylor-Hood pair of meshes is needed
    # This mesh will be for the velocities
    vmesh = Mesh()
    insert!(vmesh, T6block(A, A, N, N), "velocity")
    # This mesh will be used for the pressures
    pmesh = Mesh()
    insert!(pmesh, T6toT3(baseincrel(vmesh, "velocity")), "pressure")
  return vmesh, pmesh
end

function assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
    function integrateK!(ass, geom, elit, qpit, mu)
        for el in elit
            uxel = iterator(el, 1)
            uyel = iterator(el, 2)
            pel = iterator(el, 3)
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                for j in 1:nedof
                    DBj = D * B(gradN[j], c[j])
                    for i in 1:nedof
                        Bi = B(gradN[i], c[i])
                        v = dot(DBj, Bi) * JxW
                        asstolma!(el, i, j, v)
                    end
                end
            end
            assemble!(ass, lma(el)...)
        end
        return ass
    end

    elit = AggregateFEIterator([FEIterator(uxfesp), FEIterator(uyfesp), FEIterator(pfesp)])
    qpit = QPIterator(fesp, (kind = :default, npts = 3,))
    geom = geometry(fesp.mesh)
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    @time integrateK!(ass, geom, elit, qpit, D)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    @time KT = K * U
    @time U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function run()
    vmesh, pmesh = genmesh()
    # Velocity spaces
    uxfesp = FESpace(Float64, vmesh, FEH1_T6(), 1)
    uyfesp = FESpace(Float64, vmesh, FEH1_T6(), 1)
    locs = geometry(vmesh)
    inflate = A / N / 100
    # Part of the boundary that is immovable
    boxes = [[0.0 A 0.0 0.0], [0.0 0.0 0.0 A], [A A 0.0 A]]
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            setebc!(uxfesp, 0, i, 1, 0.0)
            setebc!(uyfesp, 0, i, 1, 0.0)
        end
    end
    # The lid
    uxbar = 1.0
    box = [0.0 A A A]
    vl = vselect(locs; box = box, inflate = inflate)
    for i in vl
        setebc!(uxfesp, 0, i, 1, uxbar)
        setebc!(uyfesp, 0, i, 1, 0.0)
    end
    # Pressure space
    pfesp = FESpace(Float64, pmesh, FEH1_T6(), 1)
    # Number the degrees of freedom
    numberfreedofs!(uxfesp)
    numberdatadofs!(uxfesp)
    numberfreedofs!(uyfesp, nunknowns(uxfesp)+1)
    numberdatadofs!(uyfesp)
    numberfreedofs!(pfesp, nunknowns(uyfesp)+1)
    numberdatadofs!(pfesp)
    @show ndofs(uxfesp)
    @show ndofs(uyfesp)
    @show ndofs(pfesp)
    @show nunknowns(uxfesp)
    @show nunknowns(uyfesp)
    @show nunknowns(pfesp)
    tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
    K = assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
    U = fill(0.0, ndofs(fesp))
    gathersysvec!(U, fesp)
    F = fill(0.0, tndof)
    solve!(U, K, F, tnunk)
    scattersysvec!(fesp, U)
    makeattribute(fesp, "U", 1:2)
    vtkwrite("stokes_driven_tht6-U", baseincrel(mesh), [(name = "U", allxyz = true)])
end

end

stokes_driven_tht6.run()
