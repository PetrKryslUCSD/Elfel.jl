"""
    th_p2_p1_veclap_alt

The manufactured-solution colliding flow example from Elman et al 2014. The
Taylor-Hood formulation with quadratic triangles for the velocity and continuous
pressure on linear triangles.

This implementation is an alternative: for the velocity, a single finite element
space with multiple components (2) is used instead of multiple finite element
spaces.

The formulation is the one derived in Donea, Huerta, Introduction to the finite element
method, 1993. Page 486 ff, and Elman, et al., Finite elements and fast
iterative solvers, p. 132. In other words, it is the vector Laplacian version.
"""
module th_p2_p1_veclap_alt

using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities, identty, attribute, VecAttrib
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
    # Taylor-Hood pair of meshes is needed
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

function assembleK(ufesp, pfesp, tndof, mu)
    function integrate!(ass, elits, qpits, mu)
        unedof, pnedof = ndofsperel.(elits)
        uedofcomp = edofcompnt(ufesp)
        kuu = LocalMatrixAssembler(unedof, unedof, 0.0)
        kup = LocalMatrixAssembler(unedof, pnedof, 0.0)
        for el in zip(elits...)
            uel, pel = el
            init!(kuu, eldofs(uel), eldofs(uel))
            init!(kup, eldofs(uel), eldofs(pel))
            for qp in zip(qpits...)
                uqp, pqp = qp
                Jac, J = jacjac(uel, uqp)
                JxW = J * weight(uqp)
                gradNu = bfungrad(uqp, Jac)
                Np = bfun(pqp)
                for j in 1:unedof, i in 1:unedof
                    if uedofcomp[i] == uedofcomp[j]
                        kuu[i, j] += (mu * JxW) * (dot(gradNu[i], gradNu[j]))
                    end
                end
                for j in 1:pnedof, i in 1:unedof
                    kup[i, j] += (-JxW * Np[j]) * gradNu[i][uedofcomp[i]]
                end
            end
            assemble!(ass, kuu)
            assemble!(ass, kup)
            assemble!(ass, transpose(kup))
        end
        return ass
    end

    elits = (FEIterator(ufesp), FEIterator(pfesp))
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(ufesp, qargs), QPIterator(pfesp, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    integrate!(ass, elits, qpits, mu)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    KT = K * U
    U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function evaluate_pressure_error(pfesp)
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

    elit = FEIterator(pfesp)
    qargs = (kind = :default, npts = 3,)
    qpit = QPIterator(pfesp, qargs)
    return integrate!(elit, qpit, truep)
end

function evaluate_velocity_error(ufesp)
    function integrate!(elit, qpit, trueux, trueuy)
        unedof = ndofsperel(elit)
        uedofcomp = edofcompnt(ufesp)
        E = 0.0
        for el in elit
            udofvals = eldofvals(el)
            for qp in qpit
                Jac, J = jacjac(el, qp)
                JxW = J * weight(qp)
                Nu = bfun(qp)
                uxt = trueux(location(el, qp)...)
                uyt = trueuy(location(el, qp)...)
                uxa = 0.0
                uya = 0.0
                for j in 1:unedof
                    (uedofcomp[j] == 1) && (uxa += (udofvals[j] * Nu[j]))
                    (uedofcomp[j] == 2) && (uya += (udofvals[j] * Nu[j]))
                end
                E += (JxW) * ((uxa - uxt)^2 + (uya - uyt)^2)
            end
        end
        return sqrt(E)
    end

    elit = FEIterator(ufesp)
    qargs = (kind = :default, npts = 3,)
    qpit = QPIterator(ufesp, qargs)
    return integrate!(elit, qpit, trueux, trueuy)
end

function run(N)
    vmesh, pmesh = genmesh(N)
    # Velocity spaces
    ufesp = FESpace(Float64, vmesh, FEH1_T6(), 2)
    locs = geometry(vmesh)
    inflate = A / N / 100
    # The entire boundary
    boxes = [[-A A -A -A], [-A -A -A A], [A A -A A], [-A A A A]]
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            setebc!(ufesp, 0, i, 1, trueux(locs[i]...))
            setebc!(ufesp, 0, i, 2, trueuy(locs[i]...))
        end
    end
    # Pressure space
    pfesp = FESpace(Float64, pmesh, FEH1_T3(), 1)
    atcenter = vselect(geometry(pmesh); nearestto = [0.0, 0.0])
    setebc!(pfesp, 0, atcenter[1], 1, 0.0)
    # Number the degrees of freedom
    numberdofs!(ufesp, pfesp)
    tndof = ndofs(ufesp) + ndofs(pfesp)
    tnunk = nunknowns(ufesp) + nunknowns(pfesp)
    # Assemble the coefficient matrix
    K = assembleK(ufesp, pfesp, tndof, mu)
    # p = spy(K, canvas = DotCanvas)
    # display(p)
    # Solve the system
    U = fill(0.0, tndof)
    gathersysvec!(U, [ufesp, pfesp])
    F = fill(0.0, tndof)
    solve!(U, K, F, tnunk)
    scattersysvec!([ufesp, pfesp], U)
    # Postprocessing
    makeattribute(pfesp, "p", 1)
    makeattribute(ufesp, "ux", 1)
    makeattribute(ufesp, "uy", 2)
    ep = evaluate_pressure_error(pfesp)
    ev = evaluate_velocity_error(ufesp)
    vtkwrite("th_p2_p1_veclap_alt-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("th_p2_p1_veclap_alt-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    return (ep, ev)
end

end

let
    N = 2
    for loop in 1:5
        ep, ev = th_p2_p1_veclap_alt.run(N)
        @show N, ep, ev
        N = N * 2
    end
end

