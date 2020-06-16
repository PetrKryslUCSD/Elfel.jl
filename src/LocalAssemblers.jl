module LocalAssemblers

struct LocalMatrixAssembler{IT<:Integer, T<:Number}
    row::Vector{IT}
    col::Vector{IT}
    M::Matrix{T}
end

function LocalMatrixAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}
    return LocalMatrixAssembler(fill(zero(IT), nrow*ncol), fill(zero(IT), nrow*ncol), fill(z, nrow, ncol))
end

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

function add!(a::L, i, j, v) where {L<:LocalMatrixAssembler{IT, T}} where {IT, T} 
    a.M[i, j] += v
    return a
end

struct LocalVectorAssembler{IT<:Integer, T<:Number}
    row::Vector{IT}
    V::Vector{T}
end

function LocalVectorAssembler(nrow::IT, z::T) where {IT, T}
    return LocalVectorAssembler(fill(zero(IT), nrow), fill(z, nrow))
end

function init!(a::L, rdofs) where {L<:LocalVectorAssembler{IT, T}} where {IT, T} 
    copyto!(a.row, rdofs)
    fill!(a.V, zero(T))
    return a
end

function add!(a::L, i, v) where {L<:LocalVectorAssembler{IT, T}} where {IT, T} 
    a.V[i] += v
    return a
end

end
