module FEFields

using StaticArrays

struct FEField{N, T, S, IT}
    shapecollection::S
    dofnums::Vector{SVector{N, IT}}
    isdatum::Vector{SVector{N, Bool}}
    dofvals::Vector{SVector{N, T}}
    v::Vector{T}
    function FEField(::Val{N}, ::Type{T}, ::Type{IT}, sc::S) where {N, T, S, IT}
        nv = nshapes(sc)
        z = fill(zero(IT), N)
        @show dofnums = [SVector{N}(z) for i in 1:nv]
        z = fill(false, N)
        isdatum = [z for i in 1:nv]
        z = fill(zero(T), N)
        dofvals = [SVector{N}(z) for i in 1:nv]
        return new{N, T, S, IT}(sc, dofnums, isdatum, dofvals)
    end
end



nterms(fef::FEField) = nshapes(fef.shapecollection)

"""
    numberdofs!(self::AbstractField)

Number the degrees of freedom.

The free components in the field are numbered consecutively. No effort is
made to optimize the numbering in any way. If you'd like to optimize the
numbering of the degrees of freedom, use the above form that sets the
permutation of the degrees of freedom, or the permutation of the nodes.
"""
# function numberdofs!(self::FEField{S}) where {N, T}
#     num = zero(typeof(self.nunknowns))
#     nents = nentities(self)
#     ndpe = ndofsperentity(self)
#     for i in 1:nents
#         n = MVector(self.nums(i))
#         for j in 1:ndpe
#             if !self.isknown._v[i][j] # free degree of freedom
#                 num = num + 1
#                 n[j] = num
#             end
#         end
#         self.nums[i] = n
#     end
#     self.nunknowns = num
#     for i in 1:nents
#         n = MVector(self.nums(i))
#         for j in 1:ndpe
#             if self.isknown._v[i][j] # free degree of freedom
#                 num = num + 1
#                 n[j] = num
#             end
#         end
#         self.nums[i] = n
#     end
#     return  self
# end

"""
    setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}

Set the EBCs (essential boundary conditions).

- `eid`  = entity identifier,
- `comp` = integer, which  degree of freedom (component),
- `val`  = value of type T
"""
# function setebc!(self::FEField{N, T}, eid, comp, val::T) where {N, T}
#     ik = MVector(self.isknown(eid))
#     ik[comp] = true
#     self.isknown[eid] = ik
#     d = MVector(self.dofvecs(eid))
#     d[comp] = val
#     self.dofvecs[eid] = d
#     return  self
# end

"""
    gathersysvec(self::FEField{N, T}) where {N, T}

Gather values from the field for the whole system vector.
"""
# function gathersysvec(self::FEField{N, T}) where {N, T}
#     nents = nentities(self)
#     ndpe = ndofsperentity(self)
#     vec = fill(zero(T), ndofs(self))
#     for i in 1:nents
#         en = self.dofnums(i)
#         for j in 1:ndpe
#             vec[en[j]] = self.dofvecs(i)[j]
#         end
#     end
#     return vec
# end

end
