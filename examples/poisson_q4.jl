module poisson_q4

using LinearAlgebra
using BenchmarkTools
# using InteractiveUtils
# using Profile
using MeshCore: retrieve, nrelations, nentities
using MeshMaker: Q4block
using MeshKeeper: Mesh, insert!, baseincrel
using Elfel.RefShapes: RefShapeTriangle, manifdim, manifdimv
using Elfel.FElements: FEH1_Q4, refshape, Jacobian, nbasisfuns
using Elfel.FEMeshes: FEMesh, geomattr, iterate, connectivity
using Elfel.FEFields: FEField, numberdofs!, ndofs
using Elfel.IntegDomains: IntegDomain, quadrule, jac, bfundata
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!, SysvecAssembler
using Elfel.Assemblers: LocalMatrixAssembler, LocalVectorAssembler, initialize!

A = 1.0 # length of the side of the square
thermal_conductivity =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
tempf(x) = tempf.(view(x, :, 1), view(x, :, 2))
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    mesh = Mesh()
    insert!(mesh, Q4block(A, A, N, N))
    return FEMesh(mesh, FEH1_Q4(1))
end

function assembleK(idom, T, kappa)
    function integrate!(ass, geom, T, fe, conn, qrule, bd, kappa)
        vmdim = Val(manifdim(refshape(fe)))
        nbf = nbasisfuns(fe)
        la = LocalMatrixAssembler(nbf, nbf, 0.0)
        for c in conn
            initialize!(la, T.dofnums, c)
            for qp in 1:qrule.npts
                gradNparams = bd[2][qp]
                Jac = jac(geom, c, gradNparams)
                J = Jacobian(vmdim, Jac)
                JxW = J * qrule.weights[qp]
                invJac = inv(Jac)
                for i in 1:nbf
                    gradNi = gradNparams[i] * invJac
                    for j in 1:nbf
                        gradNj = gradNparams[j] * invJac
                        assemble!(la, i, j, kappa * dot(gradNi, gradNj) * JxW)
                    end
                end
            end
            assemble!(ass, la)
        end
        return ass
    end
    geom = geomattr(idom.femesh)
    bd = bfundata(idom)
    ass = SysmatAssemblerSparse(0.0)
    conn = connectivity(idom.femesh)
    qrule = quadrule(idom)
    start!(ass, ndofs(T), ndofs(T))
    @time integrate!(ass, geom, T, idom.femesh.fe, conn, qrule, bd, kappa)
    return finish!(ass)
end

function assembleF(idom, Q)
    function integrate!(ass, geom, ir, qrule, nums, bd, fe, Q)
        vmdim = Val(manifdim(refshape(fe)))
        nbf = nbasisfuns(fe)
        la = LocalVectorAssembler(nbf, 0.0)
        for el in 1:nrelations(ir)
            conn = retrieve(ir, el)
            initialize!(la, nums, conn)
            for qp in 1:qrule.npts
                Ns = bd[1][qp]
                gradNparams = bd[2][qp]
                Jac = jac(geom, conn, gradNparams)
                J = Jacobian(vmdim, Jac)
                JxW = J * qrule.weights[qp]
                for i in 1:nbf
                    assemble!(la, i, Q * JxW)
                end
            end
            assemble!(ass, la)
        end
        return ass
    end
    geom = geometry(idom.fex)
    bd = bfundata(idom)
    ass = SysvecAssembler(0.0)
    ir = baseincrel(idom.fex.mesh)
    qrule = quadrule(idom)
    start!(ass, ndofs(idom.fex))
    @time integrate!(ass, geom, ir, qrule, idom.fex.field.nums, bd, fe(idom.fex.fesp), Q)
    return finish!(ass)
end

function run()
   femesh = genmesh()
   T = FEField(Float64, femesh) 
   numberdofs!(T)
   idom = IntegDomain(femesh, (kind = :Gauss, order = 2))
   K = assembleK(idom, T, thermal_conductivity)
   F = assembleF(idom, Q)
end

end

poisson_q4.run()
poisson_q4.run()
