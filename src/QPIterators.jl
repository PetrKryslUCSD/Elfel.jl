module QPIterators

using StaticArrays
using LinearAlgebra
using ..RefShapes: manifdim, IntegRule, quadrature
using ..RefShapes: npts, param_coords, weights
using ..FElements: refshape, ndofsperel
import ..FElements: bfun, bfungradpar, jacjac

# To do: Create the data structure for an element that has multiple degrees of
# freedom per entity. Also store the basic data for the scalar finite element 
# (single degree of freedom per entity).
function __bfundata(fe, qr)
    pc = qr.param_coords
    w  =  qr.weights
    npts = qr.npts
    MDIM = manifdim(refshape(fe))
    Ns = Vector{Float64}[];
    gradNparams = Vector{LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}[];
    for j in 1:npts
        push!(Ns, bfun(fe, pc[j,:]))
        push!(gradNparams, bfungradpar(fe, pc[j,:]))
    end
    return (Ns, gradNparams)
end

mutable struct QPIterator{FET, MDIM}
    fe::FET
    _quadr::IntegRule
    _bfundata::Vector{Vector{Float64}}
    _bfungraddata::Vector{Vector{LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}}
    _pt::Int64
end

function QPIterator(fe::FET, quadraturesettings) where {FET}
    _quadr = quadrature(refshape(fe), quadraturesettings)
    _bfundata = __bfundata(fe, _quadr)
    _pt = 0
    return QPIterator{FET, manifdim(refshape(fe))}(fe, _quadr, _bfundata[1], _bfundata[2], _pt)
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
    it._pt = state
    return it
end

bfun(it::QPIterator) = it._bfundata[it._pt]
bfungradpar(it::QPIterator) = it._bfungraddata[it._pt]
weight(it::QPIterator) = it._quadr.weights[it._pt]

end
