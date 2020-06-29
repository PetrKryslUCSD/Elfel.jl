"""
    th_p2_p1

The manufactured-solution colliding flow example from Elman et al 2014. The
Taylor-Hood formulation with quadratic triangles for the velocity and continuous
pressure on linear triangles.

The formulation is the one derived in Reddy, Introduction to the finite element
method, 1993. Page 486 ff.
"""
module th_p2_p1

using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities, identty, attribute, VecAttrib
using MeshSteward: T6block, T6toT3
using MeshSteward: Mesh, insert!, baseincrel, boundary
using MeshSteward: vselect, geometry, summary, transform
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T6, FEH1_T3, refshape, Jacobian
using Elfel.FESpaces: FESpace, ndofs, setebc!, nunknowns, doftype
using Elfel.FESpaces: numberfreedofs!, numberdatadofs!, numberdofs!
using Elfel.FESpaces: scattersysvec!, makeattribute, gathersysvec!, edofcompnt
using Elfel.FESpaces: highestfreedofnum, highestdatadofnum
using Elfel.FEIterators: FEIterator, ndofsperel, elnodes, eldofs, eldofvals
using Elfel.FEIterators: jacjac
using Elfel.QPIterators: QPIterator, bfun, bfungrad, weight, location
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
    insert!(vmesh, T6block(2 * A, 2 * A, N, N), "velocity")
    ir = baseincrel(vmesh)
    transform(ir, x -> x .- A)
    # This mesh will be used for the pressures. Notice that it needs to be
    # "compatible" with the velocity mesh in the sense that they need to share
    # the nodes at the corners of the triangles.
    pmesh = Mesh()
    insert!(pmesh, T6toT3(baseincrel(vmesh, "velocity")), "pressure")
    return vmesh, pmesh
end

function assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
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

    elits = (FEIterator(uxfesp), FEIterator(uyfesp), FEIterator(pfesp))
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(uxfesp, qargs), QPIterator(uyfesp, qargs), QPIterator(pfesp, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    integrateK!(ass, elits, qpits, mu)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    KT = K * U
    U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function evaluate_pressure_error(uxfesp, uyfesp, pfesp)
    geom = geometry(pfesp.mesh)
    ir = baseincrel(pfesp.mesh)
    p = attribute(ir.right, "p")
    pt = [truep(geom[i]...) for i in 1:length(geom)] 
    return norm(p - pt) / norm(pt)
end

function evaluate_velocity_error(uxfesp, uyfesp, pfesp)
    geom = geometry(uxfesp.mesh)
    ir = baseincrel(uxfesp.mesh)
    ux = attribute(ir.right, "ux")
    uy = attribute(ir.right, "uy")
    uxa = [ux[i][1] for i in 1:length(ux)]
    uya = [uy[i][1] for i in 1:length(uy)]
    uxt = [trueux(geom[i]...) for i in 1:length(geom)] 
    uyt = [trueuy(geom[i]...) for i in 1:length(geom)] 
    return sqrt(sum(vec(uxt .- uxa).^2 + vec(uyt .- uya).^2)) / sqrt(sum((uxt).^2 + (uyt).^2))
end

function evaluate_pressure_error(pfesp)
    function integrate!(elits, qpits, truep)
        pnedof = ndofsperel.(elits)
        E = 0.0
        for el in zip(elits...)
            pel = el
            dofvals = eldofvals(pel)
            for qp in zip(qpits...)
                pqp = qp
                Jac, J = jacjac(pel, pqp)
                JxW = J * weight(pqp)
                Np = bfun(pqp)
                pt = truep(location(pqp)...)
                pa = 0.0
                for j in 1:pnedof
                    pa += (dofvals[j] * Np[j])
                end
                E += (JxW) * (pa - pt)^2
            end
        end
        return E
    end

    elits = (FEIterator(pfesp),)
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(pfesp, qargs),)
    return integrate!(elits, qpits, truep)
end

function run(N)
    vmesh, pmesh = genmesh(N)
    # Velocity spaces
    uxfesp = FESpace(Float64, vmesh, FEH1_T6(), 1)
    uyfesp = FESpace(Float64, vmesh, FEH1_T6(), 1)
    locs = geometry(vmesh)
    inflate = A / N / 100
    # The entire boundary
    boxes = [[-A A -A -A], [-A -A -A A], [A A -A A], [-A A A A]]
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            setebc!(uxfesp, 0, i, 1, trueux(locs[i]...))
            setebc!(uyfesp, 0, i, 1, trueuy(locs[i]...))
        end
    end
    # Pressure space
    pfesp = FESpace(Float64, pmesh, FEH1_T3(), 1)
    atcenter = vselect(geometry(pmesh); nearestto = [0.0, 0.0])
    setebc!(pfesp, 0, atcenter[1], 1, 0.0)
    # Number the degrees of freedom
    numberdofs!(uxfesp, uyfesp, pfesp)
    @show tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    @show tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
    # Assemble the coefficient matrix
    K = assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
    # p = spy(K, canvas = DotCanvas)
    # display(p)
    # Solve the system
    U = fill(0.0, tndof)
    gathersysvec!(U, [uxfesp, uyfesp, pfesp])
    F = fill(0.0, tndof)
    solve!(U, K, F, tnunk)
    scattersysvec!([uxfesp, uyfesp, pfesp], U)
    # Postprocessing
    makeattribute(pfesp, "p", 1)
    makeattribute(uxfesp, "ux", 1)
    makeattribute(uyfesp, "uy", 1)
    @show evaluate_pressure_error(pfesp)
    @show evaluate_pressure_error(uxfesp, uyfesp, pfesp)
    @show evaluate_velocity_error(uxfesp, uyfesp, pfesp)
    vtkwrite("th_p2_p1-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("th_p2_p1-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    geom = geometry(pfesp.mesh)
    ir = baseincrel(pfesp.mesh)
    pt = VecAttrib([truep(geom[i]...) for i in 1:length(geom)])
    ir.right.attributes["pt"] = pt
    vtkwrite("th_p2_p1-pt", baseincrel(pmesh), [(name = "pt",), ])
    true
end

end

let
    N = 4
    for loop in 1:5
        th_p2_p1.run(N)
        N = N * 2
    end
end

