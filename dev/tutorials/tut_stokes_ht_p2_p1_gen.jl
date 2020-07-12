# # Solve the Stokes equation of colliding flow

# Synopsis: Compute the solution of the Stokes equation of incompressible
# viscous flow for a manufactured problem of colliding flow. Hood-Taylor
# triangular elements are used.

# The manufactured-solution colliding flow example from Elman et al 2014. The
# Hood-Taylor formulation with quadratic triangles for the velocity and
# continuous pressure on linear triangles.

# The formulation is the general elasticity-like scheme with
# strain-rate-displacement matrices. It can be manipulated into the one derived
# in Reddy, Introduction to the finite element method, 1993. Page 486 ff.

# The complete code is in the file [`tut_stokes_ht_p2_p1_gen.jl`](tut_stokes_ht_p2_p1_gen.jl).

# The solution will be defined  within a module in order to eliminate conflicts
# with data or functions defined elsewhere.

module tut_stokes_ht_p2_p1_gen

# We'll need some functionality from linear algebra, static arrays, and the mesh
# libraries. Some plotting will be produced to visualize structure of the
# stiffness matrix. Finally we will need the `Elfel` functionality.
using LinearAlgebra
using StaticArrays
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports
using UnicodePlots

# The boundary value problem is expressed in this weak form
# ```math
#  \int_{V}\underline{\underline{\varepsilon}}(\underline{\delta v})\; 
#  \underline{\underline{D}}\; \underline{\underline{\varepsilon}}(\underline{u})\; \mathrm{d} V
# - \int_{V} \mathrm{div}(\underline{\delta v}) p\; \mathrm{d} V = 0,\quad \forall \underline{\delta v}
# ```
# ```math
#  - \int_{V} \delta q\; \mathrm{div}(\underline{u}) p\; \mathrm{d} V = 0,\quad \forall \delta q
# ```
# Here ``\underline{\delta v}`` are the test functions in the velocity space, 
# and ``\delta q`` are the pressure test functions.

function run()
    mu = 1.0 # dynamic viscosity
    # This is the material-property matrix ``\underline{\underline{D}}``:
    D = SMatrix{3, 3}(
        [2*mu  0   0
          0  2*mu  0
          0    0   mu])
    A = 1.0 # half of the length of the side of the square
    N = 100 # number of element edges per side of the square
    # These three functions define the true velocity components and the true
    # pressure.
    trueux = (x, y) -> 20 * x * y ^ 3
    trueuy = (x, y) -> 5 * x ^ 4 - 5 * y ^ 4
    truep = (x, y) -> 60 * x ^ 2 * y - 20 * y ^ 3
    
    # Construct the two meshes for the mixed method. They need to support the
    # velocity and pressure spaces.
    vmesh, pmesh = genmesh(A, N)

    # Constructive velocity space: it is a vector space with two components. The
    # degrees of freedom are real numbers (`Float64`). The velocity mesh
    # carries the finite elements of  the continuity ``H ^1``, i. e. both the
    # function values and the derivatives are square integrable. Each node
    # carries 2 degrees of freedom, hence there are two velocity components per
    # node.
    ufesp = FESpace(Float64, vmesh, FEH1_T6(), 2)

    # Now we apply the boundary conditions at the nodes around the
    # circumference. 
    locs = geometry(vmesh)
    # We use searching based on the presence of the node within a box. The
    # entire boundary will be captured within these four boxes, provided we
    # inflate those boxes with a little tolerance (we can't rely on those
    # nodes to be precisely at the coordinates given, we need to introduce some
    # tolerance).
    boxes = [[-A A -A -A], [-A -A -A A], [A A -A A], [-A A A A]]
    inflate = A / N / 100
    for box in boxes
        vl = vselect(locs; box = box, inflate = inflate)
        for i in vl
            # Remember that all  components of the velocity are known at the
            # boundary.
            setebc!(ufesp, 0, i, 1, trueux(locs[i]...))
            setebc!(ufesp, 0, i, 2, trueuy(locs[i]...))
        end
    end

    # No we construct the pressure space. It is a continuous, piecewise linear
    # space supported on a mesh of three-node triangles.
    pfesp = FESpace(Float64, pmesh, FEH1_T3(), 1)

    # The pressure in this "enclosed" flow example is only known up to a constant.
    # By setting  pressure degree of freedom at one node will make the solution
    # unique.
    atcenter = vselect(geometry(pmesh); nearestto = [0.0, 0.0])
    setebc!(pfesp, 0, atcenter[1], 1, 0.0)

    # Number the degrees of freedom. First all the free degrees of freedom are
    # numbered, both velocities and pressures. Next all the data degrees of
    # freedom are numbered, again both for the velocities and for the
    # pressures.
    numberdofs!(ufesp, pfesp)
    # The total number of degrees of freedom is now calculated.
    tndof = ndofs(ufesp) + ndofs(pfesp)
    # As is the total number of unknowns.
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
    ep = evaluate_pressure_error(pfesp, truep)
    ev = evaluate_velocity_error(ufesp, trueux, trueuy)
    vtkwrite("tut_stokes_ht_p2_p1_gen-p", baseincrel(pmesh), [(name = "p",), ])
    vtkwrite("tut_stokes_ht_p2_p1_gen-v", baseincrel(vmesh), [(name = "ux",), (name = "uy",)])
    # geom = geometry(pfesp.mesh)
    # ir = baseincrel(pfesp.mesh)
    # pt = VecAttrib([truep(geom[i]...) for i in 1:length(geom)])
    # ir.right.attributes["pt"] = pt
    # vtkwrite("tut_stokes_ht_p2_p1_gen-pt", baseincrel(pmesh), [(name = "pt",), ])
    return (ep, ev)
end

function genmesh(A, N)
    # Hood-Taylor pair of meshes is needed. The first mesh is for the
    # velocities, composed of six-node triangles. 
    vmesh = attach!(Mesh(), T6block(2 * A, 2 * A, N, N), "velocity")
    # Now translate so that the center of the square is at the origin of the
    # coordinates.
    ir = baseincrel(vmesh)
    transform(ir, x -> x .- A)
    # The second mesh is used for the pressures, and it is composed of
    # three-node triangles such that the corner nodes are shared between the
    # first and the second mesh.
    pmesh = attach!(Mesh(), T6toT3(baseincrel(vmesh, "velocity")), "pressure")
    # Return the pair of meshes
    return vmesh, pmesh
end

function assembleK(ufesp, pfesp, tndof, D)
    function integrate!(ass, elits, qpits, D)
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
    integrate!(ass, elits, qpits, D)
    return finish!(ass)
end

function solve!(U, K, F, nu)
    KT = K * U
    U[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function evaluate_pressure_error(pfesp, truep)
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

function evaluate_velocity_error(ufesp, trueux, trueuy)
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

end

using .tut_stokes_ht_p2_p1_gen
tut_stokes_ht_p2_p1_gen.run()

