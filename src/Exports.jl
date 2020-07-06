module Exports


# include("FElements.jl")
# include("FEFields.jl")
# include("FESpaces.jl")
# include("QPIterators.jl")
# include("FEIterators.jl")

###############################################################################
using ..LocalAssemblers: LocalMatrixAssembler, init!
using ..LocalAssemblers: LocalVectorAssembler, init!

export LocalMatrixAssembler, init!
export LocalVectorAssembler, init!

###############################################################################
using ..Assemblers: AbstractSysmatAssembler, SysmatAssemblerSparse
using ..Assemblers: AbstractSysvecAssembler, SysvecAssembler
using ..Assemblers: start!, assemble!, finish!

export AbstractSysmatAssembler, SysmatAssemblerSparse
export AbstractSysvecAssembler, SysvecAssembler
export start!, assemble!, finish!

###############################################################################
using ..RefShapes: AbstractRefShape, RefShapePoint, RefShapeInterval, RefShapeSquare, RefShapeCube, RefShapeTriangle, RefShapeTetrahedron
using ..RefShapes: manifdimv
using ..RefShapes: IntegRule, npts, param_coords, weights, quadrature

export AbstractRefShape, RefShapePoint, RefShapeInterval, RefShapeSquare, RefShapeCube, RefShapeTriangle, RefShapeTetrahedron, manifdimv
export IntegRule, npts, param_coords, weights, quadrature

###############################################################################
using ..FElements
using ..FElements: FE, FEData
using ..FElements: shapedesc, refshape, nfeatofdim, ndofperfeat, ndofsperel, manifdim, Jacobian, jacjac, bfun, bfungradpar
using ..FElements: FEH1_L2, FEH1_T3, FEH1_T6, FEH1_Q4, FEH1_T3_BUBBLE, FEH1_T4
using ..FElements: FEL2_Q4, FEL2_T3, FEL2_T4

export FE, FEData
export shapedesc, refshape, nfeatofdim, ndofperfeat, ndofsperel, manifdim, Jacobian, jacjac, bfun, bfungradpar
export FEH1_L2, FEH1_T3, FEH1_T6, FEH1_Q4, FEH1_T3_BUBBLE, FEH1_T4
export FEL2_Q4, FEL2_T3, FEL2_T4

###############################################################################
import ..FESpaces: FESpace
import ..FESpaces: doftype, edofmdim, edofbfnum, edofcompnt
import ..FESpaces: ndofsperel, dofnum, numberfreedofs!, numberdatadofs!, ndofs
import ..FESpaces: nunknowns, highestfreedofnum, highestdatadofnum
import ..FESpaces: numberdofs!, setebc!, gathersysvec!, scattersysvec!, makeattribute

export FESpace
export doftype, edofmdim, edofbfnum, edofcompnt
export ndofsperel, dofnum, numberfreedofs!, numberdatadofs!, ndofs
export nunknowns, highestfreedofnum, highestdatadofnum
export numberdofs!, setebc!, gathersysvec!, scattersysvec!, makeattribute

###############################################################################
using ..QPIterators: QPIterator
using ..QPIterators: bfun, bfungrad, bfungradpar, weight

export QPIterator
export bfun, bfungrad, bfungradpar, weight

###############################################################################
using ..FEIterators: FEIterator
using ..FEIterators: jacjac, location, ndofsperel, elnodes, eldofs, eldofentmdims, eldofcomps, eldofvals

export FEIterator
export jacjac, location, ndofsperel, elnodes, eldofs, eldofentmdims, eldofcomps, eldofvals

end
