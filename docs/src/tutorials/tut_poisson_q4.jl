# # Solve the heat conduction equation

# Synopsis: Compute the solution of the Poisson equation of heat conduction with a
# nonzero heat source. Quadrilateral four-node elements are used.

# The solution will be defined  within a module in order to eliminate conflicts
# with data or functions defined elsewhere.

# The complete code is in the file [`tut_poisson_q4.jl`](tut_poisson_q4.jl).

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

    # Create the finite element space to represent the temperature solution. The
    # degrees of freedom are real numbers (`Float64`), the quadrilaterals are
    # defined by the mesh, and each of the elements has the continuity ``H
    # ^1``, i. e. both the function values and the derivatives are square
    # integrable.
    Uh = FESpace(Float64, mesh, FEH1_Q4())

    # Apply the essential boundary conditions at the circumference of the square
    # domain. We find the boundary incidence relation (`boundary(mesh)`), and
    # then the list of all vertices connected by the boundary cells. The
    # function `tempf` defines the analytical temperature variation, and hence
    # for each of the vertices `i` on the boundary (they are of manifold
    # dimension  `0`), we set the component of the field (1) to the exact value
    # of the temperature at that location.
    vl = connectedv(boundary(mesh));
    locs = geometry(mesh)
    for i in vl
        setebc!(Uh, 0, i, 1, tempf(locs[i]...))
    end

    # Number the degrees of freedom, both the unknowns and the data
    # (prescribed) degrees of freedom.
    numberdofs!(Uh)
    @show ndofs(Uh), nunknowns(Uh)

    # Assemble the conductivity matrix and the vector of the heat loads. Refer
    # to the definitional this function below.
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

    # Here we associate the values of the finite element space with the entities
    # of the mesh as an attribute.
    makeattribute(Uh, "T", 1)

    # The correctness of the solution is checked by comparing the values at the
    # vertices.
    checkcorrectness(Uh, tempf)

    # The attribute can now be written out for visualization into a VTK file.
    vtkwrite("q4-T", baseincrel(mesh), [(name = "T",)])

    true # return success
end

# The domain is a square, meshed with quadrilateral elements. The function
# `Q4block` creates an incidence relation that defines the quadrilateral
# element shapes by the vertices connected into the shapes. This incidence
# relation is then attached to the mesh and the mesh is returned.
function genmesh(A, N)
    conn = Q4block(A, A, N, N)
    return attach!(Mesh(), conn)
end

# This function constructs the left-hand side coefficient matrix, conductivity
# matrix, as a sparse matrix, and a vector of the heat loads due to the
# internal heat generation rate `Q`.
function assembleKF(Uh, kappa, Q)
    # This function evaluates the integrals. The key to making this calculation
    # efficient is type stability. All the arguments coming in must have
    # concrete types. This is why this function is a subfunction: the function
    # barrier allows for all arguments to be resolved to concrete types.
    function integrate!(am, av, elit, qpit, kappa, Q)
        nedof = ndofsperel(elit)
        # The local assemblers are just like matrices or vectors
        ke = LocalMatrixAssembler(nedof, nedof, 0.0) 
        fe = LocalVectorAssembler(nedof, 0.0)
        for el in elit # Loop over all elements
            init!(ke, eldofs(el), eldofs(el)) # zero out elementwise matrix
            init!(fe, eldofs(el)) # and vector
            for qp in qpit # Now loop over the quadrature points
                Jac, J = jacjac(el, qp) # Calculate the Jacobian matrix, Jacobian
                gradN = bfungrad(qp, Jac) # Evaluate the spatial gradients
                JxW = J * weight(qp) # elementary volume
                N = bfun(qp) # Basis function values at the quadrature point
                # This double loop evaluates the element wise conductivity
                # matrix and the heat load vector precisely as the formula in
                # the weak form  dictates.
                for j in 1:nedof
                    for i in 1:nedof
                        ke[i, j] += dot(gradN[i], gradN[j]) * (kappa * JxW)
                    end
                    fe[j] += N[j] * Q * JxW
                end
            end
            # Assemble the calculated contributions from this element
            assemble!(am, ke)
            assemble!(av, fe)
        end
        return am, av # Return the updated assemblers
    end

    # First we create the element iterator. We can go through all the elements that
    # define the domain of integration using this iterator. Each time a new element
    # is accessed, some data are precomputed such as the element degrees of
    # freedom.
    elit = FEIterator(Uh)
    # This is the quadrature point iterator. We know that the elements are
    # quadrilateral, which makes the Gauss integration rule the obvious choice.
    # We also select order 2 for accuracy.
    qpit = QPIterator(Uh, (kind = :Gauss, order = 2))
    # Next we create assemblers, one for the sparse system matrix and one for
    # the system vector.
    am = start!(SysmatAssemblerSparse(0.0), ndofs(Uh), ndofs(Uh))
    av = start!(SysvecAssembler(0.0), ndofs(Uh))
    # Now we call the integration function. The assemblers are modified inside
    # this function...
    @time integrate!(am, av, elit, qpit, kappa, Q)
    # ...so that when the integration is done, we can materialize the sparse
    #    matrix and the vector and return them.
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

# The correctness can be checked in various ways. Here we calculate the mean
# deviation of the calculated temperatures at the nodes relative to the exact
# values of the temperature.
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

# The module can now be used. 
using .tut_poisson_q4
tut_poisson_q4.run()
