module FEFields

using StaticArrays
using MeshKeeper: Mesh, increl, basecode
using MeshCore: nshapes, indextype, VecAttrib, AbsAttrib

struct FEField{S}
    shapecollection::S
end

function FEField(::Type{T}, sc::S) where {T, S}
    nv = nshapes(sc)
    @show dofnums =  VecAttrib([zero(typeof(nshapes(sc))) for i in 1:nv])
    @show typeof(dofnums) <: AbsAttrib
    @show typeof(sc.attributes)
    sc.attributes["dofnums"] = dofnums
    isdatum =  VecAttrib([false for i in 1:nv])
    sc.attributes["isdatum"] = isdatum
    dofvals =  VecAttrib([zero(T) for i in 1:nv])
    sc.attributes["dofvals"] = dofvals
    return FEField(sc)
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
