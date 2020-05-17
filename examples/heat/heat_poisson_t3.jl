module heat_poisson_t3

using MeshCore: retrieve, nrelations, nentities
using MeshMaker: T3block
using MeshKeeper: Mesh, insert!, baseincrel
using Elfel.RefShapes: RefShapeTriangle, manifdim, manifdimv
using Elfel.FElements: FEH1_T3, refshape, Jacobian, nbasisfuns
using Elfel.FESpaces: FESpace, fe
using Elfel.FEExpansions: FEExpansion, numberdofs!, geometry, ndofs
using Elfel.IntegDomains: IntegDomain, quadrule, jac, bfundata
using Elfel.Assemblers: SysmatAssemblerSparse, start!, finish!, assemble!, local_assembler, fill_dofs!
using LinearAlgebra
using BenchmarkTools
using InteractiveUtils
using Profile

A = 1.0 # length of the side of the square
thermal_conductivity =  1.0; # conductivity matrix
Q = -6.0; # internal heat generation rate
tempf(x, y) =(1.0 + x^2 + 2.0 * y^2);#the exact distribution of temperature
tempf(x) = tempf.(view(x, :, 1), view(x, :, 2))
N = 1000;# number of subdivisions along the sides of the square domain

function genmesh()
    conn = T3block(A, A, N, N)
    mesh = Mesh()
    insert!(mesh, conn)
    return mesh
end

function assembleK(idom)
    function integrate!(ass, geom, ir, qrule, nums, bd, fe)
        vmdim = Val(manifdim(refshape(fe)))
        nbf = nbasisfuns(fe)
        rs, cs, vs = local_assembler(nbf, nbf, 0.0)
        for el in 1:nrelations(ir)
            conn = retrieve(ir, el)
            fill!(vs, 0.0)
            fill_dofs!(rs, cs, nums, conn)
            for qp in 1:qrule.npts
                gradNparams = bd[2][qp]
                Jac = jac(geom, conn, gradNparams)
                J = Jacobian(vmdim, Jac)
                JxW = J * qrule.weights[qp]
                invJac = inv(Jac)
                k = 1
                for i in 1:nbf
                    gradNi = gradNparams[i] * invJac
                    for j in 1:nbf
                        gradNj = gradNparams[j] * invJac
                        vs[k] += dot(gradNi, gradNj) * JxW 
                        k = k + 1
                    end
                end
            end
            assemble!(ass, rs, cs, vs)
        end
        return ass
    end
    geom = geometry(idom.fex)
    bd = bfundata(idom)
    ass = SysmatAssemblerSparse(0.0)
    ir = baseincrel(idom.fex.mesh)
    qrule = quadrule(idom)
    start!(ass, ndofs(idom.fex), ndofs(idom.fex))
    @time integrate!(ass, geom, ir, qrule, idom.fex.field.nums, bd, fe(idom.fex.fesp))
    
    return finish!(ass)
end

function run()
    mesh = genmesh()
    fesp = FESpace((FEH1_T3(1), 1))
    fex = FEExpansion(mesh, fesp)
    numberdofs!(fex)
    idom = IntegDomain(fex, (kind = :default, npts = 1))
    K = assembleK(idom)
end

end

heat_poisson_t3.run()
heat_poisson_t3.run()
