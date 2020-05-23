module FESpaces

using StaticArrays
using MeshCore
using MeshCore: nshapes, indextype, nrelations, retrieve
using MeshKeeper: Mesh, baseincrel
using ..FElements: nfeatofdim, ndofsperfeat, ndofsperelem

mutable struct FEField{N, T, IT}
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}
    nunknowns::Int64
    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, ns) where {N, T, IT}
        z = fill(zero(IT), N)
        dofnums = [SVector{N}(z) for i in 1:ns]
        z = fill(false, N)
        isdatum = [z for i in 1:ns]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:ns]
        return new{N, T, IT}(sc, dofnums, isdatum, dofvals)
    end
end

struct FESpace{FET}
    fe::FET
    mesh::Mesh
    function FESpace(fe::FET, mesh) where {FET}
        return new{FET}(fe, mesh)
    end
end

struct FEIterator{FES, IR, G}
    fesp::FES
    _bir::IR
    _geom::G
    _dofs::Vector{Int64}
    _nodes::Vector{Int64}
    
    function FEIterator(fesp::FES) where {FES}
        _bir = baseincrel(fesp.mesh)
        nd = ndofsperelem(fesp.fe)
        _geom = MeshCore.attribute(_bir.right, "geom")
        _dofs = zeros(Int64, nd)
        nn = nfeatofdim(fesp.fe, 0)
        _nodes = zeros(Int64, nn) 
        return new{FES, typeof(_bir), typeof(_geom)}(fesp, _bir, _geom, _dofs, _nodes)
    end
end

function update!(i::FEIterator, state)
    #= TODO =#
    copyto!(i._nodes, retrieve(i._bir, state))
    return i
end

function Base.iterate(i::FEIterator, state = 1)
    if state > nrelations(i._bir)
        return nothing
    else
        return (update!(i, state), state+1)
    end
end
Base.length(i::FEIterator)  = nrelations(i._bir)

function numberdofs!(self::FESpace) 
end

end
