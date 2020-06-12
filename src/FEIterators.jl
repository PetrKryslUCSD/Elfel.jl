module FEIterators

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, IncRel, shapedesc
using MeshSteward: Mesh, baseincrel, increl
using ..RefShapes: manifdim, manifdimv
using ..FElements: refshape, nfeatofdim
import ..FElements: jacjac
using ..FEFields: FEField, ndofsperterm
using ..FESpaces: FESpace, doftype
import ..FESpaces: ndofsperel
using ..QPIterators: QPIterator, bfungradpar

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
struct FEIterator{FES, IR, G, IT, T, V, IR0, IR1, IR2, IR3, F0, F1, F2, F3}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{IT}
    _entmdim::Vector{IT}
    _dofcomp::Vector{IT}
    _nodes::Vector{IT}
    _irs0::Union{Nothing, IR0}
    _irs1::Union{Nothing, IR1}
    _irs2::Union{Nothing, IR2}
    _irs3::Union{Nothing, IR3}
    _fld0::Union{Nothing, F0}
    _fld1::Union{Nothing, F1}
    _fld2::Union{Nothing, F2}
    _fld3::Union{Nothing, F3}
    _lma::_LocalMatrixAssembler{IT, T}
    _lva::_LocalVectorAssembler{IT, T}
    _manifdimv::V
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, ndofsperel(fesp))
        _entmdim = zeros(Int64, ndofsperel(fesp))
        _dofcomp = zeros(Int64, ndofsperel(fesp))
        _nodes = zeros(Int64, nfeatofdim(fesp.fe, 0)) 
        _fld0 = nothing
        _fld1 = nothing
        _fld2 = nothing
        _fld3 = nothing
        _irs0 = nothing
        _irs1 = nothing
        _irs2 = nothing
        _irs3 = nothing
        0 in keys(fesp._irsfields) && (_irs0 = fesp._irsfields[0][1]; _fld0 = fesp._irsfields[0][2])
        1 in keys(fesp._irsfields) && (_irs1 = fesp._irsfields[1][1]; _fld1 = fesp._irsfields[1][2])
        2 in keys(fesp._irsfields) && (_irs2 = fesp._irsfields[2][1]; _fld2 = fesp._irsfields[2][2])
        3 in keys(fesp._irsfields) && (_irs3 = fesp._irsfields[3][1]; _fld3 = fesp._irsfields[3][2])
        _lma = _LocalMatrixAssembler(ndofsperel(fesp), zero(doftype(fesp)))
        _lva = _LocalVectorAssembler(ndofsperel(fesp), zero(doftype(fesp)))
        _manifdimv = Val(manifdim(refshape(fesp.fe)))
        p = 1
        _irs0 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs0, _fld0))
        _irs1 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs1, _fld1))
        _irs2 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs2, _fld2))
        _irs3 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs3, _fld3))
        return new{FES, typeof(_bir), typeof(_geom), eltype(_dofs), doftype(fesp), typeof(_manifdimv), typeof(_irs0), typeof(_irs1), typeof(_irs2), typeof(_irs3), typeof(_fld0), typeof(_fld1), typeof(_fld2), typeof(_fld3)}(fesp, _bir, _geom, _dofs, _entmdim, _dofcomp, _nodes, _irs0, _irs1, _irs2, _irs3, _fld0, _fld1, _fld2, _fld3, _lma, _lva, _manifdimv)
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

"""
    ndofsperel(it::FEIterator)

Retrieve the number of degrees of freedom per element.
"""
ndofsperel(it::FEIterator) = ndofsperel(it.fesp)

"""
    eldofs(it::FEIterator)

Retrieve the vector of the element degrees of freedom
"""
eldofs(it::FEIterator) = it._dofs

"""
    elnodes(it::FEIterator)

Retrieve the vector of the nodes of the element.
"""
elnodes(it::FEIterator) = it._nodes

"""
    eldofentmdims(it::FEIterator)

Retrieve the vector of the entity dimensions for each element degree of freedom.
"""
eldofentmdims(it::FEIterator) = it._entmdim

"""
    eldofcomps(it::FEIterator)

Retrieve the vector of the component numbers for each element degree of freedom.
"""
eldofcomps(it::FEIterator) = it._dofcomp

function _init_e_d!(mdim, dofcomp, p, ir, fl)
    ndpt = ndofsperterm(fl)
    for k in 1:nentities(ir, 1)
        for i in 1:ndpt
            mdim[p] = MeshCore.manifdim(shapedesc(ir.right))
            dofcomp[p] = i
            p = p + 1
        end
    end
    return p
end

function _storedofs!(d, p, e, ir, fl)
    ndpt = ndofsperterm(fl)
    for k in 1:nentities(ir, e)
        gk = retrieve(ir, e, k)
        for i in 1:ndpt
            d[p] = fl.dofnums[gk][i]
            p = p + 1
        end
    end
    return p
end

function _update!(it::FEIterator, state)
    copyto!(it._nodes, retrieve(it._bir, state))
    p = 1
    it._irs0 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs0, it._fld0))
    it._irs1 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs1, it._fld1))
    it._irs2 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs2, it._fld2))
    it._irs3 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs3, it._fld3))
    _initlma!(it)
    _initlva!(it)
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

"""
    lma(it::FEIterator)

Retrieve the local matrix assembly data.
"""
lma(it::FEIterator) = (it._lma.row, it._lma.col, it._lma.M)

"""
    asstolma!(it::FEIterator, i, j, v) 

Assemble scalar `v` into the row `i` and column `j` of the local matrix.
"""
function asstolma!(it::FEIterator, i, j, v) 
    it._lma.M[i, j] += v
    return nothing
end


function _initlva!(it::FEIterator)
    copyto!(it._lva.row, it._dofs)
    fill!(it._lva.V, zero(eltype(it._lva.V)))
    return it
end

function asstolva!(it::FEIterator, i, v) 
    it._lva.V[i] += v
    return it
end

"""
    lva(it::FEIterator)

Retrieve the local vector assembly data.
"""
lva(it::FEIterator) = (it._lva.row, it._lva.V)


"""
    jacjac(it::QPIterator)

Compute the Jacobian matrix and the Jacobian determinant.

At the current integration point.
"""
function jacjac(it::FEIterator, qpit::QPIterator)
    return jacjac(it.fesp.fe, it._geom, it._nodes, qpit._scalbfungrad_ps[qpit._pt])
end

end
