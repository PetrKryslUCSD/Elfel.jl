module IntegDomains

using StaticArrays
using Elfel.RefShapes: manifdim, IntegRule, quadrature
using Elfel.FElements: refshape, bfun, bfundpar, nbasisfuns
using Elfel.FESpaces: FESpace, multiplicity, fe
using Elfel.Fields: FEField, ndofsperentity, nentities, numberdofs!, setebc!, gathersysvec
using Elfel.FEExpansions: FEExpansion

struct IntegDomain{FET, GT}
    fex::FEExpansion{FET, GT}
    _quadr::IntegRule
    _bfundata
end

function __bfundata(fex, qr)
    pc = qr.param_coords
    w  =  qr.weights
    npts = qr.npts
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    FET = fe(fex.fesp)
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

function IntegDomain(fex::FEExpansion{FET, GT}, quadraturesettings) where {FET, GT}
	_quadr = quadrature(refshape(fe(fex.fesp)), quadraturesettings)
    _bfundata = __bfundata(fex, _quadr)
    return IntegDomain(fex, _quadr, _bfundata)
end

quadrule(idom::IntegDomain) = idom._quadr

bfundata(idom::IntegDomain) = idom._bfundata


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
