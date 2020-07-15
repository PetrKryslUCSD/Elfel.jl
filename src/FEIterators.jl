module FEIterators

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, nentities, retrieve, IncRel, shapedesc
using MeshSteward: Mesh, baseincrel, increl
using ..RefShapes: manifdim, manifdimv
using ..FElements: refshape, nfeatofdim
import ..FElements: jacjac
using ..FESpaces.FEFields: FEField, ndofsperterm
using ..FESpaces: FESpace, doftype
import ..FESpaces: ndofsperel
using ..QPIterators: QPIterator, bfungradpar

#= TODO is it more natural to have access to the geometry from the font element space or from the iterator? =#

"""
    FEIterator{FES, IR, G, IT, T, V, IR0, IR1, IR2, IR3, F0, F1, F2, F3}

Type of finite element iterator. Parameterized with the types of
- `FES`: finite element space,
- `IR`: base incidence relation of the mesh, 
- `G`: type of the geometry attribute, 
- `IT`: type of integer indices, such as the  numbers of nodes and degrees of freedom, 
- `T`: type of the degree of freedom value (real double, complex float, ... ), 
- `V`: `Val` representation of the manifold dimension of the base relation elements, 
- `IR0`, `IR1`, `IR2`, `IR3`: types of incidence relations with which degrees
  of freedom are associated in the finite element space, for each of the
  manifolds dimensions 0, 1, 2, 3, 
- `F0`, `F1`, `F2`, `F3`: types of fields with which degrees
  of freedom are associated in the finite element space, for each of the
  manifolds dimensions 0, 1, 2, 3.
"""
struct FEIterator{FES, IR, G, IT, T, V, IR0, IR1, IR2, IR3, F0, F1, F2, F3}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{IT}
    _dofvals::Vector{T}
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
    _manifdimv::V
    _state::Vector{IT}
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, ndofsperel(fesp))
        _dofvals = zeros(doftype(fesp), ndofsperel(fesp))
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
        _manifdimv = Val(manifdim(refshape(fesp.fe)))
        p = 1
        _irs0 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs0, _fld0))
        _irs1 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs1, _fld1))
        _irs2 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs2, _fld2))
        _irs3 != nothing && (p = _init_e_d!(_entmdim, _dofcomp, p, _irs3, _fld3))
        return new{FES, typeof(_bir), typeof(_geom), eltype(_dofs), doftype(fesp), typeof(_manifdimv), typeof(_irs0), typeof(_irs1), typeof(_irs2), typeof(_irs3), typeof(_fld0), typeof(_fld1), typeof(_fld2), typeof(_fld3)}(fesp, _bir, _geom, _dofs, _dofvals, _entmdim, _dofcomp, _nodes, _irs0, _irs1, _irs2, _irs3, _fld0, _fld1, _fld2, _fld3, _manifdimv, [0])
    end
end

"""
    Base.iterate(it::FEIterator, state = 1)

Advance the iterator to the next entity.

The nodes of the finite element are cached, as is a vector of all the degrees
of freedom represented on the element.
"""
function Base.iterate(it::FEIterator, state = 1)
    if state > nrelations(it._bir)
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end

"""
    Base.length(it::FEIterator) 

Number of elements represented by this iterator.
"""
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

Each degree of freedom is associated with some entity of the finite element:
vertices, edges, faces, and so on. This vector records the dimension of the
manifold entity with which each degree of freedom is associated.
"""
eldofentmdims(it::FEIterator) = it._entmdim

"""
    eldofcomps(it::FEIterator)

Retrieve the vector of the component numbers for each element degree of freedom.

If multiple copies of the finite element are referenced in the finite element
space, each copy is referred to as component.
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

function _storedofvals!(d, p, e, ir, fl)
    ndpt = ndofsperterm(fl)
    for k in 1:nentities(ir, e)
        gk = retrieve(ir, e, k)
        for i in 1:ndpt
            d[p] = fl.dofvals[gk][i]
            p = p + 1
        end
    end
    return p
end

function _update!(it::FEIterator, state)
    it._state[1] = state
    copyto!(it._nodes, retrieve(it._bir, state))
    p = 1
    it._irs0 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs0, it._fld0))
    it._irs1 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs1, it._fld1))
    it._irs2 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs2, it._fld2))
    it._irs3 != nothing && (p = _storedofs!(it._dofs, p, state, it._irs3, it._fld3))
    return it
end

"""
    eldofvals(it::FEIterator)

Provide access to vector of element degrees of freedom.
"""
function eldofvals(it::FEIterator)
    p = 1
    it._irs0 != nothing && (p = _storedofvals!(it._dofvals, p, it._state[1], it._irs0, it._fld0))
    it._irs1 != nothing && (p = _storedofvals!(it._dofvals, p, it._state[1], it._irs1, it._fld1))
    it._irs2 != nothing && (p = _storedofvals!(it._dofvals, p, it._state[1], it._irs2, it._fld2))
    it._irs3 != nothing && (p = _storedofvals!(it._dofvals, p, it._state[1], it._irs3, it._fld3))
    return it._dofvals
end

"""
    jacjac(it::FEIterator, qpit::QPIterator)

Compute the Jacobian matrix and the Jacobian determinant.

The finite element iterator cooperates with the quadrature point iterator here
to compute the Jacobian at the current integration point.
"""
function jacjac(it::FEIterator, qpit::QPIterator)
    return jacjac(it.fesp.fe, it._geom, it._nodes, qpit._geomscalbfungrad_ps[qpit._pt])
end

"""
    location(it::FEIterator, qpit::QPIterator)

Calculate the location of the quadrature point.
"""
function location(it::FEIterator, qpit::QPIterator)
    # @show it._nodes
    # @show qpit._geomscalbfuns[qpit._pt]
    n = it._nodes[1]
    loc = it._geom[n] * qpit._geomscalbfuns[qpit._pt][1]
    for i in 2:length(it._nodes)
        n = it._nodes[i]
        loc += it._geom[n] * qpit._geomscalbfuns[qpit._pt][i]
    end
    return loc
end

end
