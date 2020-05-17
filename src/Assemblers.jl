module Assemblers

using SparseArrays: sparse
using LinearAlgebra: diag

"""
    LocalAssembler{IT<:Integer, T<:Number}

Type of local assembler for a square matrix.
"""
mutable struct LocalAssembler{IT<:Integer, T<:Number}
    row::Vector{IT}
    col::Vector{IT}
    val::Vector{T}
    cursor::IT
end

function LocalAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}
    return LocalAssembler(fill(zero(IT), nrow*ncol), fill(zero(IT), nrow*ncol), fill(zero(T), nrow*ncol), zero(IT))
end

"""
    initialize!(locass, nums, conn)

Initialize a local assembler.

This needs to be done once for each finite element.
"""
function initialize!(locass, nums, conn)
    nbf = length(conn)
    k = 1
    for j in 1:nbf
        gj = nums(conn[j])[1]
        for i in 1:nbf
            gi = nums(conn[i])[1]
            locass.row[k] = gi
            locass.col[k] = gj
            k = k + 1
        end
    end
    locass.cursor = 1
    fill!(locass.val, zero(eltype(locass.val)))
    return locass
end

"""
    updatev!(locass, v)

Update the current value in the local assembler.

The cursor in the local assembler is moved forward.
"""
function updatev!(locass, v)
    locass.val[locass.cursor] += v
    locass.cursor = locass.cursor + 1
    return locass
end

"""
    nextqp!(locass)

Reset the cursor to 1 for the next quadrature point.
"""
nextqp!(locass) = (locass.cursor = 1)

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
    assemble!(self::SysmatAssemblerSparse{T}, rs::AbstractVector{IT}, cs::AbstractVector{IT}, vs::AbstractVector{T}) where {IT<:Integer, T<:Number}

Assemble the triple of the row numbers, column numbers, and values.

See the local assembler.
"""
function assemble!(self::SysmatAssemblerSparse{T}, rs::AbstractVector{IT}, cs::AbstractVector{IT}, vs::AbstractVector{T}) where {IT<:Integer, T<:Number}
    append!(self.row, rs)
    append!(self.col, cs)
    append!(self.val, vs)
    return self
end

"""
    assemble!(self::SysmatAssemblerSparse{T}, gi, ke::AbstractMatrix{T}) where {T}

Assemble a square symmetric matrix.
"""
function assemble!(self::SysmatAssemblerSparse{T}, gi::AbstractVector{IT}, ke::AbstractMatrix{T}) where {IT, T}
    n = length(gi)
    append!(self.val, ke)
    @inbounds for j in 1:n
        append!(self.row, gi)
        for i in 1:n
            push!(self.col, gi[j])
        end
    end
end

function assemble!(self::SysmatAssemblerSparse{T}, la::LocalAssembler{IT, T}) where {IT<:Integer, T<:Number}
    append!(self.row, la.row)
    append!(self.col, la.col)
    append!(self.val, la.val)
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

Assemble an elementwise vector.

The method assembles a column element vector using the vector of degree of
freedom numbers for the rows.
"""
function assemble!(self::SysvecAssembler{T}, i, val::T) where {T<:Number}
    self.val[i] = self.val[i] + val;
end

"""
    finish!(self::SysvecAssembler)

Make the global vector.
"""
function finish!(self::SysvecAssembler)
  return deepcopy(self.val);
end


function local_assembler(nrow, ncol, z::T) where {T}
	return fill(0, nrow*ncol), fill(0, nrow*ncol), fill(0.0, nrow*ncol)
end

function fill_dofs!(rs, cs, nums, conn)
	nbf = length(conn)
	k = 1
	for j in 1:nbf
		gj = nums(conn[j])[1]
		for i in 1:nbf
			gi = nums(conn[i])[1]
			rs[k] = gi
			cs[k] = gj
			k = k + 1
		end
	end
	return rs, cs
end


end
