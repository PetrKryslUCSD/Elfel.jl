module LocalAssemblers



struct LocalMatrixAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 2}
    row::Matrix{IT}
    col::Matrix{IT}
    M::Matrix{T}
end

function LocalMatrixAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}
    return LocalMatrixAssembler(fill(zero(IT), nrow, ncol), fill(zero(IT), nrow, ncol), fill(z, nrow, ncol))
end

Base.IndexStyle(::Type{<:LocalMatrixAssembler}) = IndexLinear()
Base.size(a::A) where {A<:LocalMatrixAssembler} =  size(a.M)
Base.getindex(a::A, i::Int, j::Int) where {A<:LocalMatrixAssembler} = a.M[i, j]
Base.setindex!(a::A, v, i::Int, j::Int) where {A<:LocalMatrixAssembler} =  (a.M[i, j] = v)  

function init!(a::L, rdofs, cdofs) where {L<:LocalMatrixAssembler{IT, T}} where {IT, T}
    nr = length(rdofs)
    nc = length(cdofs)
    k = 1
    for j in 1:nc
        gj = cdofs[j]
        for i in 1:nr
            gi = rdofs[i]
            a.row[k] = gi
            a.col[k] = gj
            k += 1
        end
    end
    fill!(a.M, zero(T))
    return a
end

struct LocalVectorAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 1}
    row::Vector{IT}
    V::Vector{T}
end

function LocalVectorAssembler(nrow::IT, z::T) where {IT, T}
    return LocalVectorAssembler(fill(zero(IT), nrow), fill(z, nrow))
end

Base.IndexStyle(::Type{<:LocalVectorAssembler}) = IndexLinear()
Base.size(a::A) where {A<:LocalVectorAssembler} =  size(a.V)
Base.getindex(a::A, i::Int) where {A<:LocalVectorAssembler} = a.V[i]
Base.setindex!(a::A, v, i::Int) where {A<:LocalVectorAssembler} =  (a.V[i] = v)  

function init!(a::L, rdofs) where {L<:LocalVectorAssembler{IT, T}} where {IT, T} 
    copyto!(a.row, rdofs)
    fill!(a.V, zero(T))
    return a
end

end
