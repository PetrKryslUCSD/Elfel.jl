"""
    q1_q0

The manufactured-solution colliding flow example from Elman et al 2014. The
linear quadrilaterals for the velocity and discontinuous
pressure.

The formulation is the one derived in Reddy, Introduction to the finite element
method, 1993. Page 486 ff.
"""
module q1_q0

using LinearAlgebra
using StaticArrays
using MeshCore
using MeshCore: retrieve, nrelations, nentities, ir_identity, attribute
using MeshSteward: Q4block
using MeshSteward: Mesh, attach!, baseincrel, boundary, increl
using MeshSteward: vselect, geometry, summary, transform
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEL2_Q4, FEH1_Q4, refshape, Jacobian
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

mu = 1.0 # dynamic viscosity
A = 1.0 # length of the side of the square
trueux = (x, y) -> 20 * x * y ^ 3
trueuy = (x, y) -> 5 * x ^ 4 - 5 * y ^ 4
truep = (x, y) -> 60 * x ^ 2 * y - 20 * y ^ 3

function genmesh(N)
    # This mesh will be both for the velocities and for the pressure
    mesh = Mesh()
    attach!(mesh, Q4block(2 * A, 2 * A, N, N), "velocity+pressure")
    ir = baseincrel(mesh)
    transform(ir, x -> x .- A)
    eidir = ir_identity(ir)
    attach!(mesh, eidir)
    return mesh
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
                Jac, J = jacjac(uxel, uxqp) # L2 pressure: must integrate over the velocity mesh
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
    qargs = (kind = :default, order = 2,)
    qpits = (QPIterator(uxfesp, qargs), QPIterator(uyfesp, qargs), QPIterator(pfesp, qargs))
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, tndof, tndof)
    @time integrateK!(ass, elits, qpits, mu)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    @time KT = K * U
    @time U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

centroid(p, c) = begin
    (p[c[1]] + p[c[2]] + p[c[3]] + p[c[4]]) / 4
end

function evaluate_error(uxfesp, uyfesp, pfesp)
    geom = geometry(pfesp.mesh)
    bir = baseincrel(pfesp.mesh)
    points = MeshCore.attribute(bir.right, "geom")
    centroids = [centroid(points, retrieve(bir, i)) for i in 1:nrelations(bir)] 
    p = attribute(increl(pfesp.mesh, (2, 2)).right, "p")
    pt = [truep(centroids[i]...) for i in 1:length(centroids)] 
    return norm(p - pt) / norm(pt)
end

function run(N)
    mesh = genmesh(N)
    # Velocity spaces
    uxfesp = FESpace(Float64, mesh, FEH1_Q4(), 1)
    uyfesp = FESpace(Float64, mesh, FEH1_Q4(), 1)
    locs = geometry(mesh)
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
    pfesp = FESpace(Float64, mesh, FEL2_Q4(), 1)
    atcenter = vselect(geometry(mesh); nearestto = [0.0, 0.0])
    setebc!(pfesp, 2, atcenter[1], 1, 0.0)
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
    @show evaluate_error(uxfesp, uyfesp, pfesp)
    vtkwrite("q1_q0-p", baseincrel(mesh), [(name = "p",), ])
    vtkwrite("q1_q0-v", baseincrel(mesh), [(name = "ux",), (name = "uy",)])
    true
end

end


let
    N = 4
    for loop in 1:5
        q1_q0.run(N)
        N = N * 2
    end
end

