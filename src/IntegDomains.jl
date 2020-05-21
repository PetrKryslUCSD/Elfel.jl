module IntegDomains

using StaticArrays
using Elfel.RefShapes: manifdim, IntegRule, quadrature
using Elfel.FElements: refshape, bfun, bfundpar, nbasisfuns
using Elfel.FEMeshes: FEMesh

"""
    IntegDomain{FET, GT}

Integration domain.
"""
struct IntegDomain
    femesh::FEMesh
    _quadr::IntegRule
    _bfundata
end

function __bfundata(fe, qr)
    pc = qr.param_coords
    w  =  qr.weights
    npts = qr.npts
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    FET = fe
    NBFPE = nbasisfuns(FET)
    MDIM = manifdim(refshape(FET))
    Ns = SVector{NBFPE}[];
    TMP = SVector{MDIM}(zeros(MDIM))'
    gradNparams = Vector{typeof(TMP)}[];
    gradNs = Vector{typeof(TMP)}[];
    for j in 1:npts
        push!(Ns, bfun(FET, pc[j,:]))
        push!(gradNparams, bfundpar(FET, pc[j,:]))
        push!(gradNs, bfundpar(FET, pc[j,:]))
    end
    return (Ns, gradNparams, gradNs)
end

"""
    IntegDomain(fex::FEExpansion{FET, GT}, quadraturesettings) where {FET, GT}

Create integration domain from  a finite element expansion and numerical
quadrature settings.
"""
function IntegDomain(femesh::FEMesh, quadraturesettings) 
	_quadr = quadrature(refshape(femesh.fe), quadraturesettings)
    _bfundata = __bfundata(femesh.fe, _quadr)
    return IntegDomain(femesh, _quadr, _bfundata)
end

"""
    quadrule(idom::IntegDomain)

Return the quadrature rule.
"""
quadrule(idom::IntegDomain) = idom._quadr

"""
    bfundata(idom::IntegDomain)

Return the basis function data.
"""
bfundata(idom::IntegDomain) = idom._bfundata

"""
    jac(locs, conn, gradNpar)

Compute the Jacobian matrix.
"""
function jac(locs, conn, gradNpar)
    NBFPE = length(gradNpar)
    j = 1
    J = locs[conn[j]] * gradNpar[j]
    @inbounds for j in 2:NBFPE
        J += locs[conn[j]] * gradNpar[j]
    end
    return J
end

end # module
