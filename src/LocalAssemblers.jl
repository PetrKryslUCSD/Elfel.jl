module LocalAssemblers


"""
    LocalMatrixAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 2}

Type of "local" matrix assembler.

Local is to be understood in the sense of in the context of a single finite
element. So a local matrix is an elementwise matrix which is computed entry by
entry. Then it can be assembled into the global matrix in one shot.
"""
struct LocalMatrixAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 2}
    row::Matrix{IT}
    col::Matrix{IT}
    M::Matrix{T}
end

"""
    LocalMatrixAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}

Create a local matrix assembler, given the number of rows and columns, and the
value to which the matrix should be initialized.
"""
function LocalMatrixAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}
    return LocalMatrixAssembler(fill(zero(IT), nrow, ncol), fill(zero(IT), nrow, ncol), fill(z, nrow, ncol))
end

"""
    Base.IndexStyle(::Type{<:LocalMatrixAssembler})

The data storage is assumed to be consumed best by one-dimensional traversal
through the columns in a linear fashion.
"""
Base.IndexStyle(::Type{<:LocalMatrixAssembler}) = IndexLinear()

"""
    Base.size(a::A) where {A<:LocalMatrixAssembler}

The size is the tuple of number of rows and number of columns.
"""
Base.size(a::A) where {A<:LocalMatrixAssembler} =  size(a.M)

"""
    Base.getindex(a::A, i::Int, j::Int)

Only access to a single entry of the matrix is  provided.
"""
Base.getindex(a::A, i::Int, j::Int) where {A<:LocalMatrixAssembler} = a.M[i, j]

"""
    Base.setindex!(a::A, v, i::Int, j::Int) where {A<:LocalMatrixAssembler}

Only access to a single entry of the matrix is  provided.
"""
Base.setindex!(a::A, v, i::Int, j::Int) where {A<:LocalMatrixAssembler} =  (a.M[i, j] = v)  

"""
    init!(a::L, rdofs, cdofs) where {L<:LocalMatrixAssembler{IT, T}} where {IT, T}

Initialize the  local assembler with the global degrees of freedom in the rows and columns.

The two arrays, `rdofs`, `cdofs`, define the global degree of freedom numbers
for the element. The data matrix is zeroed out. 

This function needs to be called for each new finite element.
"""
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

"""
    LocalVectorAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 1}


Type of "local" vector assembler.

Local is to be understood in the sense of in the context of a single finite
element. So a local vector is an elementwise vector which is computed entry by
entry. Then it can be assembled into the global vector in one shot.
"""
struct LocalVectorAssembler{IT<:Integer, T<:Number} <: AbstractArray{T, 1}
    row::Vector{IT}
    V::Vector{T}
end

"""
    LocalVectorAssembler(nrow::IT, z::T) where {IT, T}

Create a local vector assembler, given the number of entries, and the
value to which the vector should be initialized.
"""
function LocalVectorAssembler(nrow::IT, z::T) where {IT, T}
    return LocalVectorAssembler(fill(zero(IT), nrow), fill(z, nrow))
end

"""
    Base.IndexStyle(::Type{<:LocalVectorAssembler})

Only linear access is provided.
"""
Base.IndexStyle(::Type{<:LocalVectorAssembler}) = IndexLinear()

"""
    Base.size(a::A) where {A<:LocalVectorAssembler}

The size is the number of rows (entries).
"""
Base.size(a::A) where {A<:LocalVectorAssembler} =  size(a.V)

"""
    Base.getindex(a::A, i::Int) where {A<:LocalVectorAssembler}

Access is provided to a single entry of the vector.
"""
Base.getindex(a::A, i::Int) where {A<:LocalVectorAssembler} = a.V[i]

"""
    Base.setindex!(a::A, v, i::Int) where {A<:LocalVectorAssembler}

Access is provided to a single entry of the vector.
"""
Base.setindex!(a::A, v, i::Int) where {A<:LocalVectorAssembler} =  (a.V[i] = v)  

"""
    init!(a::L, rdofs) where {L<:LocalVectorAssembler{IT, T}} where {IT, T} 

Initialize the  local assembler with the global degrees of freedom in the rows.

The array `rdofs` defines the global degree of freedom numbers
for the element. The data vector is zeroed out. 

This function needs to be called for each new finite element.
"""
function init!(a::L, rdofs) where {L<:LocalVectorAssembler{IT, T}} where {IT, T} 
    copyto!(a.row, rdofs)
    fill!(a.V, zero(T))
    return a
end

end
