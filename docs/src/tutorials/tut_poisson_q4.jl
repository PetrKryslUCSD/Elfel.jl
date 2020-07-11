# # Solve the heat conduction equation

# Synopsis: Compute the solution of the Poisson equation of heat conduction with a
# nonzero heat source. Quadrilateral four-node elements are used.

# The solution will be defined  within a module in order to eliminate conflicts
# with data or functions defined elsewhere.

# The complete code is in the file [tut_poisson_q4.jl](tut_poisson_q4.jl).

module tut_poisson_q4

# We'll need some functionality from linear algebra, and the mesh libraries.

using LinearAlgebra
using MeshCore.Exports
using MeshSteward.Exports
using Elfel.Exports

# This is the top level function. 
function run()
    # Input parameters:
    A = 1.0 # length of the side of the square
    kappa =  1.0; # thermal conductivity of the material
    Q = -6.0; # internal heat generation rate
    tempf(x, y) =(1.0 + x^2 + 2.0 * y^2); # the exact distribution of temperature
    N = 1000; # number of element edges along the sides of the square domain

    # Generate the computational mesh.
    mesh = genmesh(A, N)

    # Create the finite element space to represent the temperature solution.
    Uh = FESpace(Float64, mesh, FEH1_Q4())

    # Apply the essential boundary conditions at the circumference of the square
    # domain.
    bir = boundary(mesh);
    vl = connectedv(bir);
    locs = geometry(mesh)
    for i in vl
        setebc!(Uh, 0, i, 1, tempf(locs[i]...))
    end

    # Number the degrees of freedom, both the unknowns and the data
    # (prescribed) degrees of freedom.
    numberdofs!(Uh)
    @show ndofs(Uh), nunknowns(Uh)

    # Assemble the conductivity matrix and the vector of the heat loads.
    K, F = assembleKF(Uh, kappa, Q)

    # This is a vector to hold all degrees of freedom in the system.
    T = fill(0.0, ndofs(Uh))
    # Here we collect the data degrees of freedom (the known values).
    gathersysvec!(T, Uh)

    # The system of linear algebraic equations is solved.
    solve!(T, K, F, nunknowns(Uh))

    # The values of all the degrees of freedom can now be introduced into the
    # finite element space.
    scattersysvec!(Uh, T)

    # Here we associate the values of the finite element with the entities of
    # the mesh as an attribute.
    makeattribute(Uh, "T", 1)

    # The correctness of the solution is checked by comparing the values at the
    # vertices.
    checkcorrectness(Uh, tempf)

    # The attribute can now be written out for visualization into a VTK file.
    vtkwrite("q4-T", baseincrel(mesh), [(name = "T",)])

    true # return success
end

# The domain is a square, meshed with quadrilateral elements.
function genmesh(A, N)
    conn = Q4block(A, A, N, N)
    return attach!(Mesh(), conn)
end

function assembleKF(Uh, kappa, Q)
    function integrate!(am, av, elit, qpit, kappa, Q)
        nedof = ndofsperel(elit)
        ke = LocalMatrixAssembler(nedof, nedof, 0.0)
        fe = LocalVectorAssembler(nedof, 0.0)
        for el in elit
            init!(ke, eldofs(el), eldofs(el))
            init!(fe, eldofs(el))
            for qp in qpit
                Jac, J = jacjac(el, qp)
                gradN = bfungrad(qp, Jac)
                JxW = J * weight(qp)
                N = bfun(qp)
                for j in 1:nedof
                    for i in 1:nedof
                        ke[i, j] += dot(gradN[i], gradN[j]) * (kappa * JxW)
                    end
                    fe[j] += N[j] * Q * JxW
                end
            end
            assemble!(am, ke)
            assemble!(av, fe)
        end
        return am, av
    end

    elit = FEIterator(Uh)
    qpit = QPIterator(Uh, (kind = :Gauss, order = 2))
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Uh), ndofs(Uh))
    av = start!(SysvecAssembler(0.0), ndofs(Uh))

    @time integrate!(am, av, elit, qpit, kappa, Q)

    return finish!(am), finish!(av)
end

# The linear algebraic system is solved by partitioning. The vector `T` is
# initially all zero, except in the degrees of freedom which are prescribed as
# nonzero. Therefore the product of the conductivity matrix and the vector `T`
# are the heat loads due to nonzero essential boundary conditions. To this we
# add the vector of heat loads due to the internal heat generation rate. The
# submatrix of the heat conduction matrix corresponding to the free degrees of
# freedom (unknowns), `K[1:nu, 1:nu]` is then used to solve for the unknowns `T
# [1:nu]`.
function solve!(T, K, F, nu)
    @time KT = K * T
    @time T[1:nu] = K[1:nu, 1:nu] \ (F[1:nu] - KT[1:nu])
end

function checkcorrectness(Uh, tempf)
    geom = geometry(Uh.mesh)
    ir = baseincrel(Uh.mesh)
    T = attribute(ir.right, "T")
    std = 0.0
    for i in 1:length(T)
        std += abs(T[i][1] - tempf(geom[i]...))
    end
    @show (std / length(T)) <= 1.0e-9
end

end # module

using .tut_poisson_q4
tut_poisson_q4.run()
