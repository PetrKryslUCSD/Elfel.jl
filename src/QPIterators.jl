module QPIterators

using StaticArrays
using LinearAlgebra
using ..RefShapes: manifdim, IntegRule, quadrature
using ..RefShapes: npts, param_coords, weights
using ..FElements: refshape, ndofsperel
import ..FElements: bfun, bfungradpar, jacjac
using ..FESpaces: edofbfnum

# To do: Create the data structure for an element that has multiple degrees of
# freedom per entity. Also store the basic data for the scalar finite element 
# (single degree of freedom per entity).
function __bfundata(fesp, qr)
    bfnum = edofbfnum(fesp)
    nedof = length(bfnum)
    pc = qr.param_coords
    w  =  qr.weights
    npts = qr.npts
    MDIM = manifdim(refshape(fesp.fe))
    # First construct vectors of vectors of basis function values and matrices
    # of base function gradients in parametric coordinates for the scalar
    # element
    scalNs = Vector{Float64}[];
    scalgradNps = Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}[];
    for j in 1:npts
        push!(scalNs, bfun(fesp.fe, pc[j,:]))
        push!(scalgradNps, bfungradpar(fesp.fe, pc[j,:]))
    end
    scalgradN = similar(scalgradNps[1])
    # Now extend it to the vector element
    Ns = Vector{Float64}[];
    gradNps = Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}[];
    for j in 1:npts
        N = [scalNs[j][bfnum[k]] for k in 1:nedof]
        gN = [scalgradNps[j][bfnum[k]] for k in 1:nedof]
        push!(Ns, N)
        push!(gradNps, gN)
    end
    gradN = similar(gradNps[1])
    return (bfnum, scalNs, scalgradNps, scalgradN, Ns, gradNps, gradN)
end

mutable struct QPIterator{FES, MDIM}
    fesp::FES
    _quadr::IntegRule
    _bfnum::Vector{Int64}
    _scalbfuns::Vector{Vector{Float64}}
    _scalbfungrad_ps::Vector{Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}}
    _scalbfungrads::Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}
    _bfuns::Vector{Vector{Float64}}
    _bfungrad_ps::Vector{Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}}
    _bfungrads::Vector{Adjoint{Float64,SArray{Tuple{MDIM},Float64,1,MDIM}}}
    _pt::Int64
end

function QPIterator(fesp::FES, quadraturesettings) where {FES}
    _quadr = quadrature(refshape(fesp.fe), quadraturesettings)
    bfnum, scalNs, scalgradNps, scalgradN, Ns, gradNps, gradN = __bfundata(fesp, _quadr)
    _pt = 0
    return QPIterator{FES, manifdim(refshape(fesp.fe))}(fesp, _quadr, bfnum, scalNs, scalgradNps, scalgradN, Ns, gradNps, gradN, _pt)
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

function bfun(it::QPIterator)
    return it._bfuns[it._pt]
end

function bfungradpar(it::QPIterator)
    return it._bfungrad_ps[it._pt]
end

function bfungrad(it::QPIterator, Jac)
    for j in 1:length(it._scalbfungrads)
        it._scalbfungrads[j] = it._scalbfungrad_ps[it._pt][j] / Jac
    end
    for k in 1:length(it._bfungrads)
        it._bfungrads[k] = it._scalbfungrads[it._bfnum[k]]
    end
    return it._bfungrads
end

weight(it::QPIterator) = it._quadr.weights[it._pt]

end
