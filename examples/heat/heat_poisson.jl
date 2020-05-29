module heat_poisson

import MeshCore: retrieve, nrelations, nentities
import MeshMaker: T3block
import MeshSteward: Mesh, insert!, baseincrel
import Elfel.RefShapes: RefShapeTriangle, manifdim
import Elfel.FElements: FEH1_T3, refshape, Jacobian, nbasisfuns
import Elfel.FESpaces: FESpace, fe
import Elfel.FEExpansions: FEExpansion, numberdofs!, geomval, ndofs
import Elfel.IntegDomains: IntegDomain, quadrule, jac, bfundata
import Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!
using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
# using Profile

A = 1.0 # length of the side of the square
thermal_conductivity =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = T3block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

# 1.1
function integrate!(ass, geom, ir, qrule, nums, bd, fe)
    vmdim = Val(manifdim(refshape(fe)))
    nbf = nbasisfuns(fe)
    is = fill(0, nbf*nbf)
    js = fill(0, nbf*nbf)
    vs = fill(0.0, nbf*nbf)
    for el in 1:nrelations(ir)
        conn = retrieve(ir, el)
        for qp in 1:qrule.npts
            gradNparams = bd[2][qp]
            Jac = jac(geom, conn, gradNparams)
            J = Jacobian(vmdim, Jac)
            JxW = J * qrule.weights[qp]
            invJac = inv(Jac)
            k = 1
            for i in 1:nbf
                gi = nums(conn[i])[1]
                gradNi = gradNparams[i] * invJac
                for j in 1:nbf
                    gj = nums(conn[j])[1]
                    gradNj = gradNparams[j] * invJac
                    v = dot(gradNi, gradNj) * JxW 
                    is[k] = gi
                    js[k] = gj
                    vs[k] = v
                    k = k + 1
                end
            end
            append!(ass.row, is)
            append!(ass.col, js)
            append!(ass.val, vs)
        end
    end
    return ass
end

# # 1.8
# function integrate!(ass, geom, ir, qrule, nums, bd, fe)
#     vmdim = Val(manifdim(refshape(fe)))
#     nbf = nbasisfuns(fe)
#     for el in 1:nrelations(ir)
#         conn = retrieve(ir, el)
#         for qp in 1:qrule.npts
#             gradNparams = bd[2][qp]
#             Jac = jac(geom, conn, gradNparams)
#             J = Jacobian(vmdim, Jac)
#             JxW = J * qrule.weights[qp]
#             invJac = inv(Jac)
#             for i in 1:nbf
#                 gi = nums(conn[i])[1]
#                 gradNi = gradNparams[i] * invJac
#                 for j in 1:nbf
#                     gj = nums(conn[j])[1]
#                     gradNj = gradNparams[j] * invJac
#                     v = dot(gradNi, gradNj) * JxW 
#                     assemble!(ass, gi, gj, v)
#                 end
#             end
#         end
#     end
#     return ass
# end

# 1.49
# function integrate!(ass, geom, ir, qrule, nums, bd, fe)
#     vmdim = Val(manifdim(refshape(fe)))
#     nbf = nbasisfuns(fe)
#     ke = fill(0.0, nbf, nbf)
#     gi = fill(0, nbf)
#     for el in 1:nrelations(ir)
#         conn = retrieve(ir, el)
#         for qp in 1:qrule.npts
#             gradNparams = bd[2][qp]
#             Jac = jac(geom, conn, gradNparams)
#             J = Jacobian(vmdim, Jac)
#             JxW = J * qrule.weights[qp]
#             invJac = inv(Jac)
#             for i in 1:nbf
#                 gi[i] = nums(conn[i])[1]
#                 gradNi = gradNparams[i] * invJac
#                 for j in 1:nbf
#                     gradNj = gradNparams[j] * invJac
#                     ke[i, j] = dot(gradNi, gradNj) * JxW 
#                 end
#             end
#             assemble!(ass, gi, ke)
#         end
#     end
#     return ass
# end

# 1.8
# function integrate!(ass, geom, ir, qrule, nums, bd, fe)
#     vmdim = Val(manifdim(refshape(fe)))
#     nbf = nbasisfuns(fe)
#     for el in 1:nrelations(ir)
#         conn = retrieve(ir, el)
#         for qp in 1:qrule.npts
#             gradNparams = bd[2][qp]
#             Jac = jac(geom, conn, gradNparams)
#             J = Jacobian(vmdim, Jac)
#             JxW = J * qrule.weights[qp]
#             invJac = inv(Jac)
#             for i in 1:nbf
#                 gi = nums(conn[i])[1]
#                 gradNi = gradNparams[i] * invJac
#                 for j in 1:nbf
#                     gj = nums(conn[j])[1]
#                     gradNj = gradNparams[j] * invJac
#                     v = dot(gradNi, gradNj) * JxW 
#                     assemble!(ass, gi, gj, v)
#                 end
#             end
#         end
#     end
#     return ass
# end

# 1.85
# function integrate!(ass, geom, ir, qrule, nums, bd, fe)
#     vmdim = Val(manifdim(refshape(fe)))
#     for el in 1:nrelations(ir)
#         conn = retrieve(ir, el)
#         for qp in 1:qrule.npts
#             gradNparams = bd[2][qp]
#             Jac = jac(geom, conn, gradNparams)
#             J = Jacobian(vmdim, Jac)
#             JxW = J * qrule.weights[qp]
#             invJac = inv(Jac)
#             gradNs = bd[3][qp]
#             for k in 1:nbasisfuns(fe)
#                 gradNs[k] = gradNparams[k] * invJac
#             end
#             # gradN = [gradNparams[i] * invJac for i in 1:nbasisfuns(fe)] 
#             for i in 1:nbasisfuns(fe)
#                 gi = nums(conn[i])[1]
#                 # gradNi = gradNparams[i] * invJac
#                 for j in 1:nbasisfuns(fe)
#                     gj = nums(conn[j])[1]
#                     # gradNj = gradNparams[j] * invJac
#                     v = dot(gradNs[i], gradNs[j]) * JxW 
#                     assemble!(ass, gi, gj, v)
#                 end
#             end
#         end
#     end
#     return ass
# end

# using StaticArrays
# using BenchmarkTools
# g = SVector{2}(rand(2))'
# J = SMatrix{2, 2}(rand(2, 2))
# @btime $g / $J # 7.299 ns (0 allocations: 0 bytes)    

function assembleK(idom)
    geom = vlocs(idom.fex)
    bd = bfundata(idom)
    ass = SysmatAssemblerSparse(0.0)
    ir = baseincrel(idom.fex.mesh)
    qrule = quadrule(idom)
    start!(ass, ndofs(idom.fex), ndofs(idom.fex))
    integrate!(ass, geom, ir, qrule, idom.fex.field.nums, bd, fe(idom.fex.fesp))

    return finish!(ass)
end

function run()
    mesh = genmesh()
    fesp = FESpace((FEH1_T3(1), 1))
    fex = FEExpansion(mesh, fesp)
    numberdofs!(fex)
    idom = IntegDomain(fex, (:default,))
    K = assembleK(idom)
end

end

heat_poisson.run()
