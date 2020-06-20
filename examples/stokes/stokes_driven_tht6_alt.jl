"""
    stokes_driven_tht6

The famous driven-cavity benchmark is solved here with Taylor-Hood combination
of quadratic and linear triangles.

The formulation is the one derived in Reddy, Introduction to the finite element
method, 1993. Page 486 ff.
"""
module stokes_driven_tht6

using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities
using MeshSteward: T6block, T6toT3
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: vselect, geometry, summary
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T6, FEH1_T3, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!, numberdofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!, edofcompnt
using Elfel.FESpaces: highestfreedofnum, highestdatadofnum
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs
using Elfel.FEIterators: jacjac
using Elfel.QPIterators: QPIterator, bfun, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!
using UnicodePlots

mu = 0.25 # dynamic viscosity
A = 1.0 # length of the side of the square
N = 100;# number of subdivisions along the sides of the square domain

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
    function integrateK!(ass, elit, qpit, mu)
        uxnedof, uynedof, pnedof = ndofsperel.(elit.is)
        kuxux = LocalMatrixAssembler(uxnedof, uxnedof, 0.0)
        kuyuy = LocalMatrixAssembler(uynedof, uynedof, 0.0)
        kuxuy = LocalMatrixAssembler(uxnedof, uynedof, 0.0)
        kuxp = LocalMatrixAssembler(uxnedof, pnedof, 0.0)
        kuyp = LocalMatrixAssembler(uynedof, pnedof, 0.0)
        for el in elit
            uxel = el[1]
            uyel = el[2]
            pel = el[3]
            init!(kuxux, eldofs(uxel), eldofs(uxel))
            init!(kuyuy, eldofs(uyel), eldofs(uyel))
            init!(kuxuy, eldofs(uxel), eldofs(uyel))
            init!(kuxp, eldofs(uxel), eldofs(pel))
            init!(kuyp, eldofs(uyel), eldofs(pel))
            for qp in qpit
                uxqp = qp[1]
                uyqp = qp[2]
                pqp = qp[3]
                Jac, J = jacjac(pel, pqp)
                JxW = J * weight(pqp)
                gradNp = bfungrad(pqp, Jac)
                gradNux = bfungrad(uxqp, Jac)
                gradNuy = bfungrad(uyqp, Jac)
                Np = bfun(pqp)
                for j in 1:uxnedof
                    for i in 1:uxnedof
                        kuxux[i, j] += (mu * JxW) * (2 * gradNux[i][1] * gradNux[j][1] + gradNux[i][2] * gradNux[j][2])
                    end
                end
                for j in 1:uynedof
                    for i in 1:uynedof
                        kuyuy[i, j] += (mu * JxW) * (gradNuy[i][1] * gradNuy[j][1] + 2 * gradNuy[i][2] * gradNuy[j][2])
                    end
                end
                for j in 1:uynedof
                    for i in 1:uxnedof
                        kuxuy[i, j] += (mu * JxW) * (gradNux[i][1] * gradNuy[j][2])
                    end
                end
                for j in 1:pnedof
                    for i in 1:uxnedof
                        kuxp[i, j] += (-JxW) * (gradNux[i][1] * Np[j])
                    end
                end
                for j in 1:pnedof
                    for i in 1:uynedof
                        kuyp[i, j] += (-JxW) * (gradNuy[i][2] * Np[j])
                    end
                end
            end
            assemble!(ass, kuxux)
            assemble!(ass, kuxuy)
            assemble!(ass, transpose(kuxuy))
            assemble!(ass, kuyuy)
            assemble!(ass, kuxp)
            assemble!(ass, transpose(kuxp))
            assemble!(ass, kuyp)
            assemble!(ass, transpose(kuyp))
        end
        return ass
    end

    elit = zip(FEIterator(uxfesp), FEIterator(uyfesp), FEIterator(pfesp))
    qargs = (kind = :default, npts = 3,)
    qpit = zip(QPIterator(uxfesp, qargs), QPIterator(uyfesp, qargs), QPIterator(pfesp, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    @time integrateK!(ass, elit, qpit, mu)
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
    pfesp = FESpace(Float64, pmesh, FEH1_T3(), 1)
    setebc!(pfesp, 0, 1, 1, 0.0)
    # Number the degrees of freedom
    numberdofs!(uxfesp, uyfesp, pfesp)
    # numberfreedofs!(uxfesp, 1)
    # numberfreedofs!(uyfesp, highestfreedofnum(uxfesp)+1)
    # numberfreedofs!(pfesp, highestfreedofnum(uyfesp)+1)
    # numberdatadofs!(uxfesp, highestfreedofnum(pfesp)+1)
    # numberdatadofs!(uyfesp, highestdatadofnum(uxfesp)+1)
    # numberdatadofs!(pfesp, highestdatadofnum(uyfesp)+1)

    @show tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    @show tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
    K = assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
    p = spy(K, canvas = DotCanvas)
    display(p)
    U = fill(0.0, tndof)
    gathersysvec!(U, uxfesp)
    gathersysvec!(U, uyfesp)
    gathersysvec!(U, pfesp)
    F = fill(0.0, tndof)
    solve!(U, K, F, tnunk)
    scattersysvec!(uxfesp, U)
    scattersysvec!(uyfesp, U)
    scattersysvec!(pfesp, U)
    makeattribute(pfesp, "p", 1)
    makeattribute(uxfesp, "ux", 1)
    makeattribute(uyfesp, "uy", 1)
    vtkwrite("stokes_driven_tht6-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("stokes_driven_tht6-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    # vtkwrite("stokes_driven_tht6-p", baseincrel(pmesh), [(name = "p", allxyz = true)])
end

end

stokes_driven_tht6.run()
