module m_stokes_driven_tht6_veclap
using Test
using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities
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
# using UnicodePlots

mu = 0.25 # dynamic viscosity
A = 1.0 # length of the side of the square
N = 100;# number of subdivisions along the sides of the square domain

function genmesh()
    # Taylor-Hood pair of meshes is needed
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

function assembleK(uxfesp, uyfesp, pfesp, tndof, mu)
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

function test()
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
    # The lid: driven in the X direction
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
    tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
    @test (tndof, tnunk) == (91003, 89402)                    
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
    vtkwrite("m_stokes_driven_tht6_veclap-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("m_stokes_driven_tht6_veclap-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    try rm("m_stokes_driven_tht6_veclap-p.vtu"); catch end
    try rm("m_stokes_driven_tht6_veclap-v.vtu"); catch end
    true
end

end
using .m_stokes_driven_tht6_veclap
m_stokes_driven_tht6_veclap.test()

module m_stokes_driven_t3b

using LinearAlgebra
using StaticArrays
using MeshCore: retrieve, nrelations, nentities, identty
using MeshSteward: T3block
using MeshSteward: Mesh, attach!, baseincrel, boundary
using MeshSteward: vselect, geometry, summary
using MeshSteward: vtkwrite
using Elfel.RefShapes: manifdim, manifdimv
using Elfel.FElements: FEH1_T3_BUBBLE, FEH1_T3, refshape, Jacobian
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
# using UnicodePlots

mu = 0.25 # dynamic viscosity
A = 1.0 # length of the side of the square
N = 100;# number of subdivisions along the sides of the square domain

function genmesh()
    # This mesh will be both for the velocities and for the pressure
    mesh = Mesh()
    attach!(mesh, T3block(A, A, N, N), "velocity+pressure")
    ir = baseincrel(mesh)
    eidir = identty(ir)
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

function test()
    mesh = genmesh()
    # Velocity spaces
    uxfesp = FESpace(Float64, mesh, FEH1_T3_BUBBLE(), 1)
    uyfesp = FESpace(Float64, mesh, FEH1_T3_BUBBLE(), 1)
    locs = geometry(mesh)
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
    pfesp = FESpace(Float64, mesh, FEH1_T3(), 1)
    setebc!(pfesp, 0, 1, 1, 0.0)
    # Number the degrees of freedom
    numberdofs!(uxfesp, uyfesp, pfesp)
    tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
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
    vtkwrite("m_stokes_driven_t3b-p", baseincrel(mesh), [(name = "p",), ])
    vtkwrite("m_stokes_driven_t3b-v", baseincrel(mesh), [(name = "ux",), (name = "uy",)])
    try rm("m_stokes_driven_t3b-p.vtu"); catch end
    try rm("m_stokes_driven_t3b-v.vtu"); catch end
    true
end

end
using .m_stokes_driven_t3b
m_stokes_driven_t3b.test()

module mstok1

using Test

"""
    th_p2_p1

The manufactured-solution colliding flow example from Elman et al 2014. The
Taylor-Hood formulation with quadratic triangles for the velocity and continuous
pressure on linear triangles.

The formulation is the one derived in Reddy, Introduction to the finite element
method, 1993. Page 486 ff.
"""
module th_p2_p1

using Test
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

function evaluate_velocity_error(uxfesp, uyfesp)
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

    elits = (FEIterator(uxfesp), FEIterator(uyfesp),)
    qargs = (kind = :default, npts = 3,)
    qpits = (QPIterator(uxfesp, qargs), QPIterator(uyfesp, qargs),)
    return integrate!(elits, qpits, trueux, trueuy)
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
    tndof = ndofs(uxfesp) + ndofs(uyfesp) + ndofs(pfesp)
    tnunk = nunknowns(uxfesp) + nunknowns(uyfesp) + nunknowns(pfesp)
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
    ep = evaluate_pressure_error(pfesp)
    ev = evaluate_velocity_error(uxfesp, uyfesp)
    # vtkwrite("th_p2_p1-p", baseincrel(pmesh), [(name = "p",), ])
    # vtkwrite("th_p2_p1-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    # geom = geometry(pfesp.mesh)
    # ir = baseincrel(pfesp.mesh)
    # pt = VecAttrib([truep(geom[i]...) for i in 1:length(geom)])
    # ir.right.attributes["pt"] = pt
    # vtkwrite("th_p2_p1-pt", baseincrel(pmesh), [(name = "pt",), ])
    (ep, ev)
end

end

function test()
    ref = [(3.5171450671095306, 0.2968271617227661),                                                                                        
    (0.5999467323539439, 0.03781189670123018),                                                   
    (0.12350320261417459, 0.004741849976722882)]         
    N = 4
    for loop in 1:3
        ep, ev = th_p2_p1.run(N)
        @test isapprox(vec([ep, ev]), vec([ref[loop]...]))
        N = N * 2
    end
end

end


using .mstok1
mstok1.test()


"""
    m_th_p2_p1_veclap_alt

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
module m_th_p2_p1_veclap_alt

using Test
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
    # vtkwrite("m_th_p2_p1_veclap_alt-p", baseincrel(pmesh), [(name = "p",), ])
    # vtkwrite("m_th_p2_p1_veclap_alt-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    return (ep, ev)
end

test() = begin
    ep, ev = m_th_p2_p1_veclap_alt.run(4)
    @test isapprox([ep, ev], [2.596076907594511, 0.3001331486426876])
end

end

using .m_th_p2_p1_veclap_alt
m_th_p2_p1_veclap_alt.test()


module m_th_p2_p1_veclap_alt_e

using Test
using LinearAlgebra
using StaticArrays
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports
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
    # vtkwrite("m_th_p2_p1_veclap_alt_e-p", baseincrel(pmesh), [(name = "p",), ])
    # vtkwrite("m_th_p2_p1_veclap_alt_e-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    return (ep, ev)
end

test() = begin
    ep, ev = m_th_p2_p1_veclap_alt_e.run(4)
    @test isapprox([ep, ev], [2.596076907594511, 0.3001331486426876])
end

end

using .m_th_p2_p1_veclap_alt_e
m_th_p2_p1_veclap_alt_e.test()
