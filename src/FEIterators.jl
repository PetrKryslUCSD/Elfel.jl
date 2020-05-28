module FEIterators

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, IncRel
using MeshSteward: Mesh, baseincrel, increl
using ..RefShapes: manifdim, manifdimv
using ..FElements: refshape, nfeatofdim, ndofsperfeat
import ..FElements: ndofsperelem, nfeatofdim, ndofsperfeat, Jacobian
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
struct FEIterator{FES, IR, G, IT, T, V}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{Int64}
    _nodes::Vector{Int64}
    _m::Vector{Int64}
    _irs::Vector{IncRel}
    _flds::Vector{FEField}
    _lma::_LocalMatrixAssembler{IT, T}
    _lva::_LocalVectorAssembler{IT, T}
    _manifdimv::V
    
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
        # _lma = _LocalMatrixAssembler(nd, zero(doftype(fesp)))
        _lma = _LocalMatrixAssembler(nd, zero(doftype(fesp)))
        _lva = _LocalVectorAssembler(nd, zero(doftype(fesp)))
        _manifdimv = Val(manifdim(refshape(fesp.fe)))
        return new{FES, typeof(_bir), typeof(_geom), eltype(_dofs), doftype(fesp), typeof(_manifdimv)}(fesp, _bir, _geom, _dofs, _nodes, _m, _irs, _flds, _lma, _lva, _manifdimv)
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
# geometry(it::FEIterator) = it._geom

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
lma(it::FEIterator) = (it._lma.row, it._lma.col, vec(it._lma.M))

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

function _jac(locs, conn, gradNpar)
    NBFPE = length(gradNpar)
    j = 1
    J = locs[conn[j]] * gradNpar[j]
    @inbounds for j in 2:NBFPE
        J += locs[conn[j]] * gradNpar[j]
    end
    return J
end

"""
    jac(it::FEIterator, gradNpar)

Compute the Jacobian matrix.
"""
function jac(it::FEIterator, gradNpar)
    return _jac(it._geom, it._nodes, gradNpar)
end

"""
    jacjac(it::FEIterator, gradNpar)

Compute the Jacobian matrix and the Jacobian determinant.
"""
function jacjac(it::FEIterator, gradNpar)
    Jac = _jac(it._geom, it._nodes, gradNpar)
    return (Jac, Jacobian(it._manifdimv, Jac))
end

end
