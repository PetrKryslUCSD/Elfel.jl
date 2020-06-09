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
    gradNps = Vector{LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}[];
    gradN = LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}};
    for j in 1:npts
        push!(Ns, bfun(fe, pc[j,:]))
        push!(gradNps, bfungradpar(fe, pc[j,:]))
    end
    gradN = similar(gradNps[1])
    return (Ns, gradNps, gradN)
end

mutable struct QPIterator{FES, MDIM}
    fesp::FES
    _quadr::IntegRule
    _bfuns::Vector{Vector{Float64}}
    _bfungrad_ps::Vector{Vector{LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}}
    _bfungrads::Vector{LinearAlgebra.Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}
    _pt::Int64
end

function QPIterator(fesp::FES, quadraturesettings) where {FES}
    _quadr = quadrature(refshape(fesp.fe), quadraturesettings)
    _bfuns = __bfundata(fesp.fe, _quadr)
    _pt = 0
    return QPIterator{FES, manifdim(refshape(fesp.fe))}(fesp, _quadr, _bfuns[1], _bfuns[2], _bfuns[3], _pt)
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

bfun(it::QPIterator) = it._bfuns[it._pt]
bfungradpar(it::QPIterator) = it._bfungrad_ps[it._pt]
bfungrad(it::QPIterator) = it._bfungrads
weight(it::QPIterator) = it._quadr.weights[it._pt]

function bfungrad(it::QPIterator, Jac)
    for j in 1:length(it._bfungrads)
        it._bfungrads[j] = it._bfungrad_ps[it._pt][j] / Jac
    end
    return it._bfungrads
end

end
