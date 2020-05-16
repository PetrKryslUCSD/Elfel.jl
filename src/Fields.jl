module Fields

import Base: setindex!
using StaticArrays

struct FEFieldVector{VT}
    _v::Vector{VT}
end

function (o::FEFieldVector{VT})(j::Int64) where {VT} 
    o._v[j]
end

function setindex!(o::FEFieldVector{VT}, X, inds...) where {VT}
    o._v[inds...] = X
    return o
end

"""
    FEField{N, T}

Structure to store degree of freedom numbers and values.
The vectors are addressed by the entity number.
"""
mutable struct FEField{N, T}
    dofvecs::FEFieldVector{SVector{N, T}}
    isknown::FEFieldVector{SVector{N, Bool}}
    nums::FEFieldVector{SVector{N, Int64}}
    nunknowns::Int64
end

function FEField{N, T}(nent::Int64) where {N, T}
    z = SVector{N, T}([zero(T) for i in 1:N])
    dofvecs = FEFieldVector([z for idx in 1:nent]) 
    z = SVector{N, Bool}([false for i in 1:N])
    isknown = FEFieldVector([z for idx in 1:nent]) 
    z = SVector{N, Int64}([0 for i in 1:N])
    nums = FEFieldVector([z for idx in 1:nent]) 
    return FEField{N, T}(dofvecs, isknown, nums, 0)
end

"""
    ndofsperentity(self::FEField{N, T}) where {N, T}

How many degrees of freedom per entity?
"""
ndofsperentity(self::FEField{N, T}) where {N, T} = N

"""
    nentities(self::FEField{N, T}) where {N, T} 

Number of entities in the field.
"""
nentities(self::FEField{N, T}) where {N, T} = length(self.dofvecs._v)

ndofs(self::FEField{N, T}) where {N, T} = nentities(self) * ndofsperentity(self)

"""
    numberdofs!(self::AbstractField)

Number the degrees of freedom.

The free components in the field are numbered consecutively. No effort is
made to optimize the numbering in any way. If you'd like to optimize the
numbering of the degrees of freedom, use the above form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
function numberdofs!(self::FEField{N, T}) where {N, T}
    num = zero(typeof(self.nunknowns))
    nents = nentities(self)
    ndpe = ndofsperentity(self)
    for i in 1:nents
        n = MVector(self.nums(i))
        for j in 1:ndpe
            if !self.isknown._v[i][j] # free degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.nums[i] = n
    end
    self.nunknowns = num
    for i in 1:nents
        n = MVector(self.nums(i))
        for j in 1:ndpe
            if self.isknown._v[i][j] # free degree of freedom
                num = num + 1
                n[j] = num
            end
        end
        self.nums[i] = n
    end
    return  self
end

"""
    setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}

Set the EBCs (essential boundary conditions).

- `eid`  = entity identifier,
- `comp` = integer, which  degree of freedom (component),
- `val`  = value of type T
"""
function setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}
    ik = MVector(self.isknown(eid))
    ik[comp] = true
    self.isknown[eid] = ik
    d = MVector(self.dofvecs(eid))
    d[comp] = val
    self.dofvecs[eid] = d
    return  self
end

"""
    gathersysvec(self::FEField{N, T}) where {N, T}

Gather values from the field for the whole system vector.
"""
function gathersysvec(self::FEField{N, T}) where {N, T}
    nents = nentities(self)
    ndpe = ndofsperentity(self)
    vec = fill(zero(T), ndofs(self))
    for i in 1:nents
        en = self.nums(i)
        for j in 1:ndpe
            vec[en[j]] = self.dofvecs(i)[j]
        end
    end
    return vec
end

# """
#     gathersysvec!(self::F, vec::FVec{T}) where {F<:AbstractField, T}

# Gather values from the field for the whole system vector.
# """
# function gathersysvec!(self::F, vec::FVec{T}) where {F<:AbstractField, T}
#     nents,dim = size(self.values)
#     @assert length(vec) == self.nfreedofs
#     for i = 1:nents
#         for j = 1:dim
#             en = self.dofnums[i,j]
#             if (en > 0) && (en <= self.nfreedofs)
#                 vec[en] = self.values[i,j]
#             end
#         end
#     end
#     return vec
# end

# """
#     gathervalues_asvec!(self::AbstractField, dest::AbstractArray{T, 1},
#         conn::CC) where {CC, T}

# Gather values from the field into a vector.

# The order is: for each node  in the connectivity, copy into the buffer all the
# degrees of freedom,  then the next node and so on.

# `dest` = destination buffer: overwritten  inside,  must be preallocated
# in the correct size
# """
# function gathervalues_asvec!(self::AbstractField, dest::AbstractArray{T, 1}, conn::CC) where {CC, T}
#     # The order of the loops matters, first i, then j
#     en::FInt = 1;
#     for i = 1:length(conn)
#         for j = 1:size(self.values,2)
#             dest[en] = self.values[conn[i],j];
#             en = en + 1;
#         end
#     end
#     return dest
# end

# """
#     gathervalues_asmat!(self::AbstractField, dest::AbstractArray{T, 2},
#         conn::CC) where {CC, T}

# Gather values from the field into a two-dimensional array.

# The order is: for each node  in the connectivity, copy into the corresponding
# row of the buffer all the degrees of freedom,  then the next node into the next
# row and so on.

# `dest` = destination buffer: overwritten  inside,  must be preallocated
# in the correct size
# """
# function gathervalues_asmat!(self::AbstractField, dest::AbstractArray{T, 2},    conn::CC) where {CC, T}
#     @inbounds for j in 1:size(self.values,2)
#         @inbounds for i in 1:length(conn)
#             dest[i, j] = self.values[conn[i], j];
#         end
#     end
#     return dest
# end

# """
#     gatherfixedvalues_asvec!(self::AbstractField, dest::AbstractArray{T, 1},
#         conn::CC) where {CC, T}

# Gather FIXED values from the field into a vector.

# The order is: for each node  in the connectivity, copy into the buffer all the
# fixed degrees of freedom,  then the next node and so on. If a degree of freedom
# is NOT fixed, the corresponding entry is  set to zero.

# `dest` = destination buffer: overwritten  inside,  must be preallocated
# in the correct size
# """
# function gatherfixedvalues_asvec!(self::AbstractField, dest::AbstractArray{T, 1},    conn::CC) where {CC, T}
#     # The order of the loops matters here! It must be i, j
#     en::FInt = 1;
#     for i = 1:length(conn)
#         for j = 1:size(self.fixed_values,2)
#             if self.is_fixed[conn[i],j] # free degree of freedom
#                 dest[en] = self.fixed_values[conn[i], j];
#             else
#                 dest[en] = 0.0
#             end
#             en = en + 1;
#         end
#     end
#     return dest
# end

# """
#     gatherfixedvalues_asmat!(self::AbstractField, dest::AbstractArray{T, 2},
#         conn::CC) where {CC, T}

# Gather FIXED values from the field into a two-dimensional array.

# The order is: for each node  in the connectivity, copy into the corresponding
# row of the buffer all the degrees of freedom,  then the next node into the next
# row and so on.  If a degree of freedom
# is NOT fixed, the corresponding entry is  set to zero.

# `dest` = destination buffer: overwritten  inside,  must be preallocated
# in the correct size
# """
# function gatherfixedvalues_asmat!(self::AbstractField, dest::AbstractArray{T, 2},    conn::CC) where {CC, T}
#     for j = 1:size(self.fixed_values,2)
#         for i = 1:length(conn)
#             if self.is_fixed[conn[i],j] # fixed degree of freedom
#                 dest[i, j] = self.fixed_values[conn[i], j];
#             else
#                 dest[i, j] = 0.0
#             end
#         end
#     end
#     return dest
# end

# function anyfixedvaluenz(self::AbstractField, conn::CC) where {CC}
#     for i = 1:length(conn)
#         for j = 1:size(self.fixed_values,2)
#             if self.is_fixed[conn[i],j] # free degree of freedom
#                 if  abs(self.fixed_values[conn[i], j]) > 0.0
#                     return true
#                 end
#             end
#         end
#     end
#     return false
# end

# """
#     gatherdofnums!(self::AbstractField, dest::A, conn::CC) where {A, CC}

# Gather dofnums from the field.

# The order is: for each node  in the connectivity, copy into the buffer all the degrees of
# freedom for that node,  then the next node  and so on.
# """
# function gatherdofnums!(self::AbstractField, dest::A, conn::CC) where {A, CC}
#     en::FInt = 1;
#     for i = 1:length(conn)
#         for j = 1:size(self.dofnums,2)
#             dest[en] = self.dofnums[conn[i],j];
#             en = en+1;
#         end
#     end
#     return dest
# end

# """
#     numberdofs!(self::AbstractField)

# Number the degrees of freedom.

# The free components in the field are numbered consecutively. No effort is
# made to optimize the numbering in any way. If you'd like to optimize the
# numbering of the degrees of freedom, use the above form that sets the
# permutation of the degrees of freedom, or the permutation of the nodes.
# """
# function numberdofs!(self::AbstractField)
#     fixed_dofnum::FInt = 0
#     nents,dim = size(self.values)
#     self.nfreedofs::FInt =0
#     for i=1:nents
#         for j=1:dim
#             if !self.is_fixed[i,j] # free degree of freedom
#                 self.nfreedofs = self.nfreedofs + 1
#                 self.dofnums[i,j] = self.nfreedofs
#             else # fixed degree of freedom: no equation
#                 self.dofnums[i,j] = fixed_dofnum
#             end
#         end
#     end
#     return  self
# end

# function _setebc!(self::AbstractField, fenid::FInt, is_fixed::Bool, comp::FInt, val::T) where {T<:Number}
#     self.is_fixed[fenid,comp] = is_fixed;
#     if self.is_fixed[fenid,comp]
#     	self.fixed_values[fenid,comp] = val;
#     else
#     	self.fixed_values[fenid,comp] = zero(T)
#     end
#     return  self
# end



# """
#     setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FInt,
#       val::FVec{T}) where {T<:Number}

# Set the EBCs (essential boundary conditions).

# `fenids`         - array of N node identifiers
# `is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
#              or released (false),
# `comp` = integer, which  degree of freedom (component),
# `val` = array of N values of type T

# Note:  Any call to `setebc!()` potentially changes the current assignment
# which degrees of freedom are free and which are fixed and therefore is
# presumed to invalidate the current degree-of-freedom numbering. In such a case
# this method sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FInt, val::FVec{T}) where {T<:Number}
#     @assert comp <= size(self.values,2) "Requested  nonexistent  degree of freedom"
#     @assert maximum(fenids) <= size(self.values,1) "Requested nonexistent node"
#     @assert size(fenids) == size(val) "Arrays of mismatched sizes"
#     for  j = 1:length(fenids)
#         _setebc!(self, fenids[j], is_fixed, comp, val[j])
#     end
#     self.nfreedofs = 0
#     fill!(self.dofnums, 0)
#     return  self
# end

# """
#     setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FInt, val::T = 0.0) where {T<:Number}

# Set the EBCs (essential boundary conditions).

# `fenids`         - array of N node identifiers
# `is_fixed` = scaler Boolean: are the degrees of freedom being fixed (true)
#              or released (false),
# `comp` = integer, which  degree of freedom (component),
# `val` = scalar of type T

# Note:  Any call to `setebc!()` potentially changes the current assignment
# which degrees of freedom are free and which are fixed and therefore is
# presumed to invalidate the current degree-of-freedom numbering. In such a case
# this method sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FInt, val::T = 0.0) where {T<:Number}
#     @assert (comp >= 1 && comp <= size(self.values,2)) "Requested  nonexistent  degree of freedom"
#     @assert maximum(fenids) <= size(self.values,1) "Requested nonexistent node"
#     @assert minimum(fenids) >= 1 "Requested nonexistent node"
#     for  j = 1:length(fenids)
#         _setebc!(self, fenids[j], is_fixed, comp, val)
#     end
#     self.nfreedofs = 0
#     fill!(self.dofnums, 0)
#     return  self
# end

# """
#     setebc!(self::AbstractField, fenids::FIntVec, comp::FInt,
#       val::FVec{T}) where {T<:Number}

# Set the EBCs (essential boundary conditions).

# `fenids` = array of N node identifiers
# `comp` = integer, which  degree of freedom (component),
# `val` = array of N values of type `T`

# Note:  Any call to `setebc!()` potentially changes the current assignment
# which degrees of freedom are free and which are fixed and therefore is
# presumed to invalidate the current degree-of-freedom numbering. In such a case
# this method sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec, comp::FInt, val::FVec{T}) where {T<:Number}
#     return setebc!(self, fenids, true, comp, val)
# end


# """
#     setebc!(self::AbstractField, fenids::FIntVec, comp::FInt, val::T=0.0) where {T<:Number}

# Set the EBCs (essential boundary conditions).

# `fenids` = array of N node identifiers
# `comp` = integer, which  degree of freedom (component),
# `val` = scalar of type `T`

# Note:  Any call to `setebc!()` potentially changes the current assignment
# which degrees of freedom are free and which are fixed and therefore is
# presumed to invalidate the current degree-of-freedom numbering. In such a case
# this method sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec, comp::FInt, val::T=0.0) where {T<:Number}
#     return setebc!(self, fenids, true, comp, val)
# end

# """
#     setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FIntVec, val::T=0.0) where {T<:Number}

# Set the EBCs (essential boundary conditions).

# `fenids` = array of N node identifiers
# `comp` = integer vector, which degree of freedom (component),
# `val` = scalar of type `T`, default is `0.0`

# Note:  Any call to `setebc!()` potentially changes the current assignment which
# degrees of freedom are free and which are fixed and therefore is presumed to
# invalidate the current degree-of-freedom numbering. In such a case this method
# sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec, is_fixed::Bool, comp::FIntVec, val::T=0.0) where {T<:Number}
# 	for j in comp
# 		setebc!(self, fenids, is_fixed, j, val)
# 	end
#     return self
# end

# """
#     setebc!(self::AbstractField, fenids::FIntVec)

# Set the EBCs (essential boundary conditions).

# Suppress all degrees of freedom at the given nodes.

# `fenids`         - array of N node identifiers

# Note:  Any call to `setebc!()` potentially changes the current assignment
# which degrees of freedom are free and which are fixed and therefore is
# presumed to invalidate the current degree-of-freedom numbering. In such a case
# this method sets `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenids::FIntVec)
#     zer = zero(eltype(self.fixed_values[1]))
#     for comp in 1:size(self.values, 2)
#         setebc!(self, fenids, true, comp, zer)
#     end
#     return self
# end

# """
#     setebc!(self::AbstractField, fenid::FInt)

# Set the EBCs (essential boundary conditions).

# Suppress all degrees of freedom at the given node.

# `fenid`         - One integer as a node identifier

# Note:  Any call to setebc!() potentially changes the current assignment
# which degrees of freedom are free and which are fixed
# and therefore is presumed to invalidate the
# current degree-of-freedom numbering. In such a case this method sets
# `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField, fenid::FInt)
#     return setebc!(self, [fenid])
# end

# """
#     setebc!(self::AbstractField)

# Set the EBCs (essential boundary conditions).

# All essential boundary conditions are CLEARED.

# Note:  Any call to setebc!() potentially changes the current assignment
# which degrees of freedom are free and which are fixed
# and therefore is presumed to invalidate the
# current degree-of-freedom numbering. In such a case this method sets
# `nfreedofs = 0`; and  `dofnums=0`.
# """
# function setebc!(self::AbstractField)
#     self.nfreedofs = 0
#     fill!(self.dofnums, 0)
#     fill!(self.is_fixed, false)
#     fill!(self.fixed_values, zero(eltype(self.fixed_values[1])))
#     return  self
# end

# """
#     applyebc!(self::AbstractField)

# Apply EBCs (essential boundary conditions).
# """
# function applyebc!(self::AbstractField)
#     nents,dim = size(self.values);
#     for i = 1:nents
#         for j = 1:dim
#             if self.is_fixed[i,j]
#                 self.values[i,j] = self.fixed_values[i,j];
#             end
#         end
#     end
#     return  self
# end

# """
#     scattersysvec!(self::AbstractField, vec::FVec{T}) where {T<:Number}

# Scatter values to the field from a system vector.
# """
# function scattersysvec!(self::AbstractField, vec::FVec{T}) where {T<:Number}
#     nents,dim = size(self.values);
#     for i = 1:nents
#         for j = 1:dim
#             dn = self.dofnums[i,j];
#             if (dn > 0) && (dn <= self.nfreedofs)
#                 self.values[i,j] = vec[dn];
#             end
#         end
#     end
#     return  self
# end

# """
#     incrscattersysvec!(self::AbstractField, vec::FVec{T}) where {T<:Number}

# Increment values of the field by scattering a system vector.
# """
# function incrscattersysvec!(self::AbstractField, vec::FVec{T}) where {T<:Number}
#     nents,dim = size(self.values);
#     for i = 1:nents
#         for j = 1:dim
#             dn = self.dofnums[i,j];
#             if (dn > 0) && (dn <= self.nfreedofs)
#                 self.values[i,j] += vec[dn];
#             end
#         end
#     end
#     return  self
# end

# """
#     prescribeddofs(uebc, u)

# Find which degrees of freedom are prescribed.
# `uebc` = field which defines the constraints (is the dof fixed and to which value),
# `u` = field which does not have the constraints applied, and serves as the source of equation numbers,
# `uebc` and `u` may be one and the same field.

# """
# function prescribeddofs(uebc::AbstractField, u::AbstractField)
#     dofnums = FInt[]
# 	prescribedvalues = eltype(uebc.values)[]
# 	nents, dim = size(uebc.values)
# 	@assert size(uebc.values) == size(u.values)
# 	for i in 1:nents
# 		for j in 1:dim
# 			if uebc.is_fixed[i,j]
# 			    dn = u.dofnums[i,j];
#                 push!(prescribedvalues, uebc.fixed_values[i,j]);
# 			    push!(dofnums, dn)
# 			end
# 		end
# 	end
# 	return dofnums, prescribedvalues
# end

end
