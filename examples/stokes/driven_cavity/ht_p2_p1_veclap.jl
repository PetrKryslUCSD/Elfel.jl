"""
    ht_p2_p1_veclap

The famous driven-cavity benchmark is solved here with Hood-Taylor combination
of quadratic and linear triangles.

The formulation is the one derived in Elman, et al., Finite elements and fast
iterative solvers, p. 132. In other words, it is the vector Laplacian version.
"""
module ht_p2_p1_veclap

using LinearAlgebra
using StaticArrays
using MeshCore: nrelations, nentities
using MeshSteward: T6block, T6toT3
using MeshSteward: Mesh, attach!, baseincrel, boundary
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
    # Hood-Taylor pair of meshes is needed
    # This mesh will be for the velocities
    vmesh = Mesh()
    attach!(vmesh, T6block(A, A, N, N), "velocity")
    # This mesh will be used for the pressures. Notice that it needs to be
        # "compatible" with the velocity mesh in the sense that they need to share
        # the nodes at the corners of the triangles.
        pmesh = Mesh()
    attach!(pmesh, T6toT3(baseincrel(vmesh, "velocity")), "pressure")
    return vmesh, pmesh
end

function assembleK(Uxh, Uyh, Ph, tndof, mu)
    function integrateK!(ass, elits, qpits, mu)
        uxnedof, uynedof, pnedof = ndofsperel.(elits)
        kuxux = LocalMatrixAssembler(uxnedof, uxnedof, 0.0)
        kuyuy = LocalMatrixAssembler(uynedof, uynedof, 0.0)
        kuxp = LocalMatrixAssembler(uxnedof, pnedof, 0.0)
        kuyp = LocalMatrixAssembler(uynedof, pnedof, 0.0)
        for el in zip(elits...)
            uxel, uyel, pel = el
            init!(kuxux, eldofs(uxel), eldofs(uxel))
            init!(kuyuy, eldofs(uyel), eldofs(uyel))
            init!(kuxp, eldofs(uxel), eldofs(pel))
            init!(kuyp, eldofs(uyel), eldofs(pel))
            for qp in zip(qpits...)
                uxqp, uyqp, pqp = qp
                Jac, J = jacjac(pel, pqp)
                JxW = J * weight(pqp)
                gradNp = bfungrad(pqp, Jac)
                gradNux = bfungrad(uxqp, Jac)
                gradNuy = bfungrad(uyqp, Jac)
                Np = bfun(pqp)
                for j in 1:uxnedof, i in 1:uxnedof
                    kuxux[i, j] += (mu * JxW) * (2 * gradNux[i][1] * gradNux[j][1] + gradNux[i][2] * gradNux[j][2])
                end
                for j in 1:uynedof, i in 1:uynedof
                    kuyuy[i, j] += (mu * JxW) * (gradNuy[i][1] * gradNuy[j][1] + 2 * gradNuy[i][2] * gradNuy[j][2])
                end
                for j in 1:pnedof, i in 1:uxnedof
                    kuxp[i, j] += (-JxW) * (gradNux[i][1] * Np[j])
                end
                for j in 1:pnedof, i in 1:uynedof
                    kuyp[i, j] += (-JxW) * (gradNuy[i][2] * Np[j])
                end
            end
            assemble!(ass, kuxux)
            assemble!(ass, kuyuy)
            assemble!(ass, kuxp)
            assemble!(ass, transpose(kuxp))
            assemble!(ass, kuyp)
            assemble!(ass, transpose(kuyp))
        end
        return ass
    end

    elits = (FEIterator(Uxh), FEIterator(Uyh), FEIterator(Ph))
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(Uxh, qargs), QPIterator(Uyh, qargs), QPIterator(Ph, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    @time integrateK!(ass, elits, qpits, mu)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    @time KT = K * U
    @time U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function run()
    vmesh, pmesh = genmesh()
    # Velocity spaces
    Uxh = FESpace(Float64, vmesh, FEH1_T6(), 1)
    Uyh = FESpace(Float64, vmesh, FEH1_T6(), 1)
    locs = geometry(vmesh)
    inflate = A / N / 100
    # Part of the boundary that is immovable
    boxes = [[0.0 A 0.0 0.0], [0.0 0.0 0.0 A], [A A 0.0 A]]
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            setebc!(Uxh, 0, i, 1, 0.0)
            setebc!(Uyh, 0, i, 1, 0.0)
        end
    end
    # The lid: driven in the X direction
    uxbar = 1.0
    box = [0.0 A A A]
    vl = vselect(locs; box = box, inflate = inflate)
    for i in vl
        setebc!(Uxh, 0, i, 1, uxbar)
        setebc!(Uyh, 0, i, 1, 0.0)
    end
    # Pressure space
    Ph = FESpace(Float64, pmesh, FEH1_T3(), 1)
    setebc!(Ph, 0, 1, 1, 0.0)
    # Number the degrees of freedom
    numberdofs!([Uxh, Uyh, Ph])
    @show tndof = ndofs(Uxh) + ndofs(Uyh) + ndofs(Ph)
    @show tnunk = nunknowns(Uxh) + nunknowns(Uyh) + nunknowns(Ph)
    # Assemble the coefficient matrix
    K = assembleK(Uxh, Uyh, Ph, tndof, mu)
    p = spy(K, canvas = DotCanvas)
    display(p)
    # Solve the system
    U = fill(0.0, tndof)
    gathersysvec!(U, [Uxh, Uyh, Ph])
    F = fill(0.0, tndof)
    solve!(U, K, F, tnunk)
    scattersysvec!([Uxh, Uyh, Ph], U)
    # Postprocessing
    makeattribute(Ph, "p", 1)
    makeattribute(Uxh, "ux", 1)
    makeattribute(Uyh, "uy", 1)
    vtkwrite("ht_p2_p1_veclap-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("ht_p2_p1_veclap-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
end

end

ht_p2_p1_veclap.run()
