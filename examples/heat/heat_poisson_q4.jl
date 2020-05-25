module heat_poisson_q4

using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
# using Profile
using MeshCore: retrieve, nrelations, nentities
using MeshMaker: Q4block
using MeshKeeper: Mesh, insert!, baseincrel
using Elfel.RefShapes: RefShapeTriangle, manifdim, manifdimv
using Elfel.FElements: FEH1_Q4, refshape, Jacobian, jac
using Elfel.FESpaces: FESpace, ndofs, numberdofs!, setebc!, nunknowns, doftype
using Elfel.FEIterators: FEIterator, asstolma!, lma, geometry, ndofsperelem, elemnodes, elemdofs
using Elfel.QPIterators: QPIterator, bfun, bfungradpar, weight
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!

A = 1.0 # length of the side of the square
thermal_conductivity =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
tempf(x) = tempf.(view(x, :, 1), view(x, :, 2))
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = Q4block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function integrateK!(ass, geom, elit, qpit)
    vmdim = Val(manifdim(refshape(elit.fesp.fe)))
    nedof = ndofsperelem(elit)
    is = fill(0, nedof*nedof)
    js = fill(0, nedof*nedof)
    vs = fill(0.0, nedof*nedof)
    for el in elit
        for qp in qpit
            gradNparams = bfungradpar(qp)
            Jac = jac(geom, elemnodes(el), gradNparams)
            J = Jacobian(vmdim, Jac)
            JxW = J * weight(qp)
            invJac = inv(Jac)
            k = 1
            for j in 1:nedof
                gradNj = gradNparams[j] * invJac
                for i in 1:nedof
                    gradNi = gradNparams[i] * invJac
                    v = dot(gradNi, gradNj) * JxW 
                    asstolma!(el, i, j, v)
                    # el._lma.M[i, j] += v
                    # is[k] = elemdofs(el)[i]
                    # js[k] = elemdofs(el)[j]
                    # vs[k] = v
                    k = k + 1
                end
            end
        end
        assemble!(ass, lma(el)...)
        # assemble!(ass, el._lma.row, el._lma.col, vec(el._lma.M))
        # assemble!(ass, is, js, vs)
    end
    return ass
end

# function f(qpit) 
#     # @show typeof(qpit)
#     gradNparams = bfungradpar(qpit)
#     # gradNparams = weight(qpit)
# end

# function f(d, j)
#     # gradNparams = bfungradpar(qpit)
#     gradNparams = d[j]
# end
# function f() 
#     # @show typeof(qpit)
#     qpit = QPIterator(fesp.fe, (kind = :Gauss, order = 2))
#     gradNparams = bfungradpar(qpit)
#     # gradNparams = weight(qpit, 1)
# end

function assembleK(fesp)
    elit = FEIterator(fesp)
    qpit = QPIterator(fesp.fe, (kind = :Gauss, order = 2))
    geom = geometry(elit)
    ass = SysmatAssemblerSparse(0.0)
    start!(ass, ndofs(fesp), ndofs(fesp))
    # @code_warntype f(qpit)
    # @code_warntype f(qpit._bfungraddata, qpit._pt[])
    @time integrateK!(ass, geom, elit, qpit)
    return finish!(ass)
end

function run()
    mesh = genmesh()
    fesp = FESpace(Float64, FEH1_Q4(1), mesh)
    numberdofs!(fesp)
    K = assembleK(fesp)
end

end

heat_poisson_q4.run()
# heat_poisson_q4.run()
