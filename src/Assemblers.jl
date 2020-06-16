module Assemblers

using SparseArrays: sparse
using LinearAlgebra
using ..LocalAssemblers: LocalMatrixAssembler, init!

"""
    AbstractSysmatAssembler

Abstract type of system-matrix assembler.
"""
abstract type AbstractSysmatAssembler end;

"""
    SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler

Type for assembling a sparse global matrix from individual entries.
"""
mutable struct SysmatAssemblerSparse{T<:Number} <: AbstractSysmatAssembler
    row::Vector{Int64}
    col::Vector{Int64}
    val::Vector{T}
    nrow::Int64; ncol::Int64;
end

"""
    SysmatAssemblerSparse(zero::T=0.0) where {T<:Number}

Construct blank system matrix assembler. The matrix entries are of type `T`.

# Example

This is how a sparse matrix is assembled from two rectangular dense matrices.
```
a = SysmatAssemblerSparse(0.0)                                                        
start!(a, 7, 7)  
m = [0.24406   0.599773    0.833404  0.0420141                                             
0.786024  0.00206713  0.995379  0.780298                                              
0.845816  0.198459    0.355149  0.224996]     
gi = [1 7 5]             
gj = [5 2 1 4]       
for j in 1:size(m, 2), i in 1:size(m, 1)
    assemble!(a, gi[i], gj[j], m[i, j])       
end  
m = [0.146618  0.53471   0.614342    0.737833                                              
0.479719  0.41354   0.00760941  0.836455                                              
0.254868  0.476189  0.460794    0.00919633                                            
0.159064  0.261821  0.317078    0.77646                                               
0.643538  0.429817  0.59788     0.958909]                                   
gi =  [2 3 1 7 5]
gj = [6 7 3 4]   
for j in 1:size(m, 2), i in 1:size(m, 1)
    assemble!(a, gi[i], gj[j], m[i, j])       
end                               
A = finish!(a) 
```
"""
function SysmatAssemblerSparse(zero::T=0.0) where {T<:Number}
    return SysmatAssemblerSparse{T}(fill(0, 0), fill(0, 0), fill(zero, 0), 0, 0)
end

"""
    start!(self::SysmatAssemblerSparse{T}, nrow, ncol) where {T<:Number}

Start the assembly of a global matrix.
"""
function start!(self::SysmatAssemblerSparse{T}, nrow, ncol) where {T<:Number}
    nexpected = Int64(round(0.00001 * max(nrow, ncol, 10000)^2))
    self.row = fill(0, 0)
    sizehint!(self.row, nexpected)
    self.col = fill(0, 0)
    sizehint!(self.col, nexpected)
    self.val = fill(zero(T), 0)
    sizehint!(self.val, nexpected)
    self.nrow = nrow;
    self.ncol = ncol;
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, r, c, v::T) where {T<:Number}

Assemble a single entry of a rectangular matrix.
"""
function assemble!(self::SysmatAssemblerSparse{T}, r, c, v::T) where {T<:Number}
    push!(self.row, r)
    push!(self.col, c)
    push!(self.val, v)
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, lma::LocalMatrixAssembler{IT, T}) where {IT<:Integer, T<:Number}

Assemble the row numbers, column numbers, and values from a local assembler.
"""
function assemble!(self::SysmatAssemblerSparse{T}, lma::LocalMatrixAssembler{IT, T}) where {IT<:Integer, T<:Number}
    append!(self.row, lma.row)
    append!(self.col, lma.col)
    append!(self.val, lma.M)
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, lma::Transpose{T,LocalMatrixAssembler{IT,T}}) where {IT<:Integer, T<:Number}

Assemble the row numbers, column numbers, and values from a local assembler.
"""
function assemble!(self::SysmatAssemblerSparse{T}, lma::Transpose{T,LocalMatrixAssembler{IT,T}}) where {IT<:Integer, T<:Number}
    append!(self.row, transpose(lma.parent.col))
    append!(self.col, transpose(lma.parent.row))
    append!(self.val, transpose(lma.parent.M))
    return self
end

"""
    finish!(self::SysmatAssemblerSparse)

Make a sparse matrix.
"""
function finish!(self::SysmatAssemblerSparse)
    return sparse(self.row, self.col, self.val, self.nrow, self.ncol);
end


"""
    AbstractSysvecAssembler

Abstract type of system vector assembler.
"""
abstract type AbstractSysvecAssembler end;

"""
    start!(self::SV,  nrow) where {SV<:AbstractSysvecAssembler, T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`nrow`= Total number of degrees of freedom.
"""
function start!(self::SV,  nrow) where {SV<:AbstractSysvecAssembler, T<:Number}
end

"""
    assemble!(self::SV, i, val::T) where {SV<:AbstractSysvecAssembler, T<:Number}

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SV, i, val::T) where {SV<:AbstractSysvecAssembler, T<:Number}
end

"""
    makevector!(self::SysvecAssembler)

Make the global vector.
"""
function finish!(self::SV) where {SV<:AbstractSysvecAssembler}
end

"""
    SysvecAssembler

Assembler for the system vector.
"""
mutable struct SysvecAssembler{T<:Number} <: AbstractSysvecAssembler
    val::Vector{T};
    ndofs::Int64
end

"""
    SysvecAssembler(zero::T=0.0) where {T<:Number}

Construct blank system vector assembler. The vector entries are of type `T`.
"""
function SysvecAssembler(zero::T=0.0) where {T<:Number}
    return SysvecAssembler([zero], 0)
end

"""
    start!(self::SysvecAssembler{T},
      nrow::FInt) where {T<:Number}

Start assembly.

The method makes the buffer for the vector assembly. It must be called before
the first call to the method assemble.

`nrow`= Total number of degrees of freedom.
"""
function start!(self::SysvecAssembler{T},  nrow::Int64) where {T<:Number}
    self.ndofs = nrow
    self.val = zeros(T, self.ndofs);
    return self
end

"""
    assemble!(self::SysvecAssembler{T}, vec::MV,
      dofnums::D) where {T<:Number, MV<:AbstractArray{T}, D<:AbstractArray{FInt}}

Assemble a single value.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssembler{T}, i, val::T) where {T<:Number}
    self.val[i] = self.val[i] + val;
    return self
end

"""
    assemble!(self::SysvecAssembler{T}, rs::AbstractVector{IT}, vs::AbstractVector{T}) where {IT<:Integer, T<:Number}

Assemble entire vector.
"""
function assemble!(self::SysvecAssembler{T}, lva) where {T<:Number}
    for i in 1:length(lva.row)
        gi = lva.row[i]
        self.val[gi] += lva.V[i];
    end
    return self
end

"""
    finish!(self::SysvecAssembler)

Make the global vector.
"""
function finish!(self::SysvecAssembler)
  return deepcopy(self.val);
end

end
