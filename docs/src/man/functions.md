# Functions

```@meta
CurrentModule = Elfel
```

## Reference shapes

```@docs
RefShapes.manifdim
RefShapes.manifdimv
RefShapes.quadrature
```


## Elements

```@docs
FElements.shapedesc
FElements.refshape
FElements.nfeatofdim
FElements.ndofperfeat
FElements.ndofsperel
FElements.manifdim
FElements.Jacobian
FElements.jacjac
FElements.bfun
FElements.bfungradpar
FElements.FEH1_L2
FElements.FEH1_T3
FElements.FEH1_Q4
FElements.FEH1_T6
FElements.FEH1_T3_BUBBLE
```

## Spaces

```@docs
FESpaces.doftype
FESpaces.edofmdim
FESpaces.edofbfnum
FESpaces.edofcompnt
FESpaces.ndofsperel(fesp::FESpaces.FESpace)
FESpaces.numberfreedofs!(fesp::FESpaces.FESpace, firstnum = 1)
FESpaces.numberdatadofs!(fesp::FESpaces.FESpace, firstnum = 0)
FESpaces.ndofs(fesp::FESpaces.FESpace)
FESpaces.nunknowns
FESpaces.highestfreedofnum(fesp::FESpaces.FESpace)
FESpaces.highestdatadofnum(fesp::FESpaces.FESpace)
FESpaces.numberdofs!
FESpaces.setebc!(fesp::FESpaces.FESpace, mid, eid, comp, val)
FESpaces.gathersysvec!(v, fesp::FESpaces.FESpace)
FESpaces.gathersysvec!(v, fesp::AbstractVector) 
FESpaces.scattersysvec!(fesp::FESpaces.FESpace, v)
FESpaces.scattersysvec!(fesp::AbstractVector, v) 
FESpaces.makeattribute
```

### Spaces/Fields

```@docs
FESpaces.FEFields.doftype
FESpaces.FEFields.dofnumtype
FESpaces.FEFields.nterms
FESpaces.FEFields.ndofsperterm
FESpaces.FEFields.ndofs
FESpaces.FEFields.setebc!
FESpaces.FEFields.numberfreedofs!
FESpaces.FEFields.numberdatadofs!
FESpaces.FEFields.freedofnums
FESpaces.FEFields.datadofnums
FESpaces.FEFields.highestfreedofnum
FESpaces.FEFields.highestdatadofnum
FESpaces.FEFields.gathersysvec!
FESpaces.FEFields.scattersysvec!
```

## Finite element iterators

```@docs
Base.iterate
FEIterators.ndofsperel(it::FEIterators.FEIterator)
FEIterators.eldofs
FEIterators.elnodes
FEIterators.eldofentmdims
FEIterators.eldofcomps
FEIterators.jacjac(it::FEIterators.FEIterator, qpit::QPIterators.QPIterator)
```

## Quadrature-point iterators

```@autodocs
Modules = [Elfel.QPIterators]
Pages   = ["QPIterators.jl"]
```

## Assemblers

```@docs
Assemblers.SysmatAssemblerSparse(zero::T=0.0) where {T<:Number}
Assemblers.start!
Assemblers.assemble!
Assemblers.finish!
Assemblers.SysvecAssembler(zero::T=0.0) where {T<:Number}
```

## Local Assemblers

```@docs
LocalAssemblers.LocalMatrixAssembler(nrow::IT, ncol::IT, z::T) where {IT, T}
LocalAssemblers.LocalVectorAssembler(nrow::IT, z::T) where {IT, T}
Base.size
Base.getindex
Base.setindex!
LocalAssemblers.init!
```

## Index

```@index
```

