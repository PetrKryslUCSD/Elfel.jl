module QPIterators

using StaticArrays
using ..RefShapes: manifdim, IntegRule, quadrature
using ..RefShapes: npts, param_coords, weights
using ..FElements: refshape, ndofsperelem
import ..FElements: bfun, bfungradpar

function __bfundata(fe, qr)
    pc = qr.param_coords
    w  =  qr.weights
    npts = qr.npts
    # Precompute basis f. values + basis f. gradients wrt parametric coor
    FET = fe
    NBFPE = ndofsperelem(FET)
    MDIM = manifdim(refshape(FET))
    Ns = SVector{NBFPE}[];
    TMP = SVector{MDIM}(zeros(MDIM))'
    gradNparams = Vector{typeof(TMP)}[];
    gradNs = Vector{typeof(TMP)}[];
    for j in 1:npts
        push!(Ns, bfun(FET, pc[j,:]))
        push!(gradNparams, bfungradpar(FET, pc[j,:]))
        push!(gradNs, bfungradpar(FET, pc[j,:]))
    end
    return (Ns, gradNparams, gradNs)
end

struct QPIterator{FET, BFT}
    fe::FET
    _quadr::IntegRule
    _bfundata::BFT
    _pt::Ref{Int64}
    
    function QPIterator(fe::FET, quadraturesettings) where {FET}
        _quadr = quadrature(refshape(fe), quadraturesettings)
        _bfundata = __bfundata(fe, _quadr)
        _pt = Ref(0)
        return new{FET, typeof(_bfundata)}(fe, _quadr, _bfundata, _pt)
    end
end

function Base.iterate(it::QPIterator, state = 1)
    if state > npts(it._quadr)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::QPIterator)  = npts(it._quadr)

function _update!(it::QPIterator, state)
    it._pt[] = state
    return it
end

bfun(it::QPIterator) = it._bfundata[1][it._pt[]]
bfungradpar(it::QPIterator) = it._bfundata[2][it._pt[]]
weight(it::QPIterator) = it._quadr.weights[it._pt[]]

end
