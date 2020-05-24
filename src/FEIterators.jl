module FEIterators

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, manifdim, IncRel
using MeshKeeper: Mesh, baseincrel, increl
using ..FElements: nfeatofdim, ndofsperfeat
import ..FElements: ndofsperelem
using ..FElements: nfeatofdim, ndofsperfeat
using ..FEFields: FEField
using ..FESpaces: FESpace, doftype

struct _LocalMatrixAssembler{IT<:Integer, T<:Number}
    row::Vector{IT}
    col::Vector{IT}
    M::Matrix{T}
end

function _LocalMatrixAssembler(nrow::IT, z::T) where {IT, T}
    ncol = nrow
    return _LocalMatrixAssembler(fill(zero(IT), nrow*ncol), fill(zero(IT), nrow*ncol), fill(z, nrow, ncol))
end

struct _LocalVectorAssembler{IT<:Integer, T<:Number}
    row::Vector{IT}
    V::Vector{T}
end

function _LocalVectorAssembler(nrow::IT, z::T) where {IT, T}
    return _LocalVectorAssembler(fill(zero(IT), nrow), fill(z, nrow))
end

#= TODO is it more natural to have access to the geometry from the font element space or from the iterator? =#
struct FEIterator{FES, IR, G}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{Int64}
    _nodes::Vector{Int64}
    _m::Vector{Int64}
    _irs::Vector{IncRel}
    _flds::Vector{FEField}
    _lma::_LocalMatrixAssembler
    _lva::_LocalVectorAssembler
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        nd = ndofsperelem(fesp.fe)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, nd)
        nn = nfeatofdim(fesp.fe, 0)
        _nodes = zeros(Int64, nn) 
        _m = Int64[]
        _irs = IncRel[]
        _flds = FEField[]
        for m in keys(fesp._irsfields)
            v = fesp._irsfields[m]
            push!(_m, m)
            push!(_irs, v[1])
            push!(_flds, v[2])
        end
        _lma = _LocalMatrixAssembler(nd, zero(doftype(fesp)))
        _lva = _LocalVectorAssembler(nd, zero(doftype(fesp)))
        return new{FES, typeof(_bir), typeof(_geom)}(fesp, _bir, _geom, _dofs, _nodes, _m, _irs, _flds, _lma, _lva)
    end
end

function Base.iterate(it::FEIterator, state = 1)
    if state > nrelations(it._bir)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::FEIterator)  = nrelations(it._bir)

ndofsperelem(it::FEIterator) = ndofsperelem(it.fesp.fe)
elemdofs(it::FEIterator) = it._dofs
elemnodes(it::FEIterator) = it._nodes

function _storedofs!(d, p, e, ir, fl)
    c = retrieve(ir, e)
    for k in 1:nentities(ir, e)
        gk = retrieve(ir, e, k)
        for i in 1:length(fl.dofnums[gk])
            d[p] = fl.dofnums[gk][i]
            p = p + 1
        end
    end
    return p
end

function _update!(it::FEIterator, state)
    copyto!(it._nodes, retrieve(it._bir, state))
    p = 1
    for i in 1:length(it._m)
        m = it._m[i]
        p = _storedofs!(it._dofs, p, state, it._irs[i], it._flds[i])
    end
    _initlma!(it)
    return it
end

function _initlma!(it)
    nd = length(it._dofs)
    k = 1
    for j in 1:nd
        gj = it._dofs[j]
        for i in 1:nd
            gi = it._dofs[i]
            it._lma.row[k] = gi
            it._lma.col[k] = gj
            k = k + 1
        end
    end
    fill!(it._lma.M, zero(eltype(it._lma.M)))
    return it
end

lma(it::FEIterator) = (it._lma.row, it._lma.col, vec(it._lma.M))

"""
    asstolma!(it::FEIterator, i, j, v) 

Assemble scalar `v` into the row `i` and column `j` of the local matrix.
"""
function asstolma!(it::FEIterator, i, j, v) 
    it._lma.M[i, j] += v
    return it
end


function _initlva!(it)
    nd = length(it._dofs)
    for i in 1:nd
        gi = it._dofs[i]
        it._lva.row[i] = gi
    end
    fill!(it._lva.V, zero(eltype(it._lva.V)))
    return it
end

function asstolva!(it::FEIterator, i, v) 
    it._lva.V[i] += v
    return it
end

end
