module FESpaces

using ..RefShapes: RefShapeInterval, RefShapeTriangle
using ..FElements: AbstractFE, refshape
using StaticArrays

"""
    FESpace

Type to represent the finite element space  for an (in general) vector degree of
freedom per entity. 
"""
struct FESpace{T<:AbstractFE}
    feinfo::Tuple{T, Int64}
end

function FESpace(feinfo::Tuple{T, Int64}) where {T}
    return FESpace(feinfo)
end

fe(self::FESpace{T}) where {T<:AbstractFE} = self.feinfo[1]
multiplicity(self::FESpace{T}) where {T<:AbstractFE} = self.feinfo[2]

end # module
