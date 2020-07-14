"""
    ht_p2_p1_gen

The manufactured-solution colliding flow example from Elman et al 2014. The
Hood-Taylor formulation with quadratic triangles for the velocity and continuous
pressure on linear triangles.

The formulation is the general elasticity-like scheme with
strain-rate-displacement matrices. It can be manipulated into the one derived in
Reddy, Introduction to the finite element method, 1993. Page 486 ff.
"""
module ht_p2_p1_gen

using LinearAlgebra
using StaticArrays
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports
using UnicodePlots

mu = 1.0 # dynamic viscosity
D = SMatrix{3, 3}(
    [2*mu 0 0
     0 2*mu 0
     0 0 mu])
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

function assembleK(ufesp, pfesp, tndof, D)
    function integrateK!(ass, elits, qpits, D)
        B = (g, k) -> (k == 1 ? 
            SVector{3}((g[1], 0, g[2])) : 
            SVector{3}((0, g[2], g[1])))
        c = edofcompnt(ufesp)
        unedof, pnedof = ndofsperel.((ufesp, pfesp))
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
                for j in 1:unedof
                    DBj = D * B(gradNu[j], c[j])
                    for i in 1:unedof
                        Bi = B(gradNu[i], c[i])
                        kuu[i, j] += dot(Bi, DBj) * (JxW)
                    end
                end
                for j in 1:pnedof, i in 1:unedof
                    kup[i, j] += (-JxW * Np[j]) * gradNu[i][c[i]]
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
    integrateK!(ass, elits, qpits, D)
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
    # Velocity space: space with two components
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
    numberdofs!([ufesp, pfesp])
    tndof = ndofs(ufesp) + ndofs(pfesp)
    tnunk = nunknowns(ufesp) + nunknowns(pfesp)
    # Assemble the coefficient matrix
    K = assembleK(ufesp, pfesp, tndof, D)
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
    vtkwrite("ht_p2_p1_gen-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("ht_p2_p1_gen-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    # geom = geometry(pfesp.mesh)
    # ir = baseincrel(pfesp.mesh)
    # pt = VecAttrib([truep(geom[i]...) for i in 1:length(geom)])
    # ir.right.attributes["pt"] = pt
    # vtkwrite("ht_p2_p1_gen-pt", baseincrel(pmesh), [(name = "pt",), ])
    return (ep, ev)
end

end

using .ht_p2_p1_gen
let
    N = 4
    for loop in 1:5
        ep, ev = ht_p2_p1_gen.run(N)
        @show N, ep, ev
        N = N * 2
    end
end

