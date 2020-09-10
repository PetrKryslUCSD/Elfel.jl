"""
    ht_p2_p1

The manufactured-solution colliding flow example from Elman et al 2014. The
Hood-Taylor formulation with quadratic triangles for the velocity and continuous
pressure on linear triangles.

The formulation is the one derived in Reddy, Introduction to the finite element
method, 1993. Page 486 ff.
"""
module ht_p2_p1

using LinearAlgebra
using StaticArrays
using MeshCore: nrelations, nentities, ir_identity, attribute, VecAttrib
using MeshSteward: T6block, T6toT3
using MeshSteward: Mesh, attach!, baseincrel, boundary
using MeshSteward: vselect, geometry, summary, transform
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T6, FEH1_T3, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!, numberdofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!, edofcompnt
using Elfel.FESpaces: highestfreedofnum, highestdatadofnum
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs, eldofvals
using Elfel.FEIterators: jacjac, location
using Elfel.QPIterators: QPIterator, bfun, bfungrad, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using Elfel.Assemblers: SysvecAssembler
using Elfel.LocalAssemblers: LocalMatrixAssembler, LocalVectorAssembler, init!
using UnicodePlots

mu = 1.0 # dynamic viscosity
A = 1.0 # half of the length of the side of the square
trueux = (x, y) -> 20 * x * y ^ 3
trueuy = (x, y) -> 5 * x ^ 4 - 5 * y ^ 4
truep = (x, y) -> 60 * x ^ 2 * y - 20 * y ^ 3

function genmesh(N)
    # Hood-Taylor pair of meshes is needed
    # This mesh will be for the velocities
    vmesh = Mesh()
    attach!(vmesh, T6block(2 * A, 2 * A, N, N), "velocity")
    ir = baseincrel(vmesh)
    transform(ir, x -> x .- A)
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
        kuxuy = LocalMatrixAssembler(uxnedof, uynedof, 0.0)
        kuxp = LocalMatrixAssembler(uxnedof, pnedof, 0.0)
        kuyp = LocalMatrixAssembler(uynedof, pnedof, 0.0)
        for el in zip(elits...)
            uxel, uyel, pel = el
            init!(kuxux, eldofs(uxel), eldofs(uxel))
            init!(kuyuy, eldofs(uyel), eldofs(uyel))
            init!(kuxuy, eldofs(uxel), eldofs(uyel))
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
                for j in 1:uynedof, i in 1:uxnedof
                    kuxuy[i, j] += (mu * JxW) * (gradNux[i][1] * gradNuy[j][2])
                end
                for j in 1:pnedof, i in 1:uxnedof
                    kuxp[i, j] += (-JxW) * (gradNux[i][1] * Np[j])
                end
                for j in 1:pnedof, i in 1:uynedof
                    kuyp[i, j] += (-JxW) * (gradNuy[i][2] * Np[j])
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

    elits = (FEIterator(Uxh), FEIterator(Uyh), FEIterator(Ph))
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(Uxh, qargs), QPIterator(Uyh, qargs), QPIterator(Ph, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    integrateK!(ass, elits, qpits, mu)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    KT = K * U
    U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function evaluate_pressure_error(Ph)
    function integrate!(elit, qpit, truep)
        pnedof = ndofsperel(elit)
        E = 0.0
        for el in elit
            dofvals = eldofvals(el)
            for qp in qpit
                Jac, J = jacjac(el, qp)
                JxW = J * weight(qp)
                Np = bfun(qp)
                pt = truep(location(el, qp)...)
                pa = 0.0
                for j in 1:pnedof
                    pa += (dofvals[j] * Np[j])
                end
                E += (JxW) * (pa - pt)^2
            end
        end
        return sqrt(E)
    end

    elit = FEIterator(Ph)
    qargs = (kind = :default, npts = 3,)
    qpit = QPIterator(Ph, qargs)
    return integrate!(elit, qpit, truep)
end

function evaluate_velocity_error(Uxh, Uyh)
    function integrate!(elits, qpits, trueux, trueuy)
        uxnedof, uynedof = ndofsperel.(elits)
        E = 0.0
        for el in zip(elits...)
            uxel, uyel = el
            uxdofvals = eldofvals(uxel)
            uydofvals = eldofvals(uyel)
            for qp in zip(qpits...)
                uxqp, uyqp = qp
                Jac, J = jacjac(uxel, uxqp)
                JxW = J * weight(uxqp)
                Np = bfun(uxqp)
                uxt = trueux(location(uxel, uxqp)...)
                uyt = trueuy(location(uyel, uyqp)...)
                uxa = 0.0
                uya = 0.0
                for j in 1:uxnedof
                    uxa += (uxdofvals[j] * Np[j])
                    uya += (uydofvals[j] * Np[j])
                end
                E += (JxW) * ((uxa - uxt)^2 + (uya - uyt)^2)
            end
        end
        return sqrt(E)
    end

    elits = (FEIterator(Uxh), FEIterator(Uyh),)
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(Uxh, qargs), QPIterator(Uyh, qargs),)
    return integrate!(elits, qpits, trueux, trueuy)
end

function run(N)
    vmesh, pmesh = genmesh(N)
    # Velocity spaces
    Uxh = FESpace(Float64, vmesh, FEH1_T6(), 1)
    Uyh = FESpace(Float64, vmesh, FEH1_T6(), 1)
    locs = geometry(vmesh)
    inflate = A / N / 100
    # The entire boundary
    boxes = [[-A A -A -A], [-A -A -A A], [A A -A A], [-A A A A]]
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            setebc!(Uxh, 0, i, 1, trueux(locs[i]...))
            setebc!(Uyh, 0, i, 1, trueuy(locs[i]...))
        end
    end
    # Pressure space
    Ph = FESpace(Float64, pmesh, FEH1_T3(), 1)
    atcenter = vselect(geometry(pmesh); nearestto = [0.0, 0.0])
    setebc!(Ph, 0, atcenter[1], 1, 0.0)
    # Number the degrees of freedom
    numberdofs!([Uxh, Uyh, Ph])
    tndof = ndofs(Uxh) + ndofs(Uyh) + ndofs(Ph)
    tnunk = nunknowns(Uxh) + nunknowns(Uyh) + nunknowns(Ph)
    # Assemble the coefficient matrix
    K = assembleK(Uxh, Uyh, Ph, tndof, mu)
    # p = spy(K, canvas = DotCanvas)
    # display(p)
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
    ep = evaluate_pressure_error(Ph)
    ev = evaluate_velocity_error(Uxh, Uyh)
    vtkwrite("ht_p2_p1-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("ht_p2_p1-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    geom = geometry(Ph.mesh)
    ir = baseincrel(Ph.mesh)
    pt = VecAttrib([truep(geom[i]...) for i in 1:length(geom)])
    ir.right.attributes["pt"] = pt
    vtkwrite("ht_p2_p1-pt", baseincrel(pmesh), [(name = "pt",), ])
    return (ep, ev)
end

end

using .ht_p2_p1
let
    N = 4
    for loop in 1:5
        ep, ev = ht_p2_p1.run(N)
        @show N, ep, ev
        N = N * 2
    end
end

