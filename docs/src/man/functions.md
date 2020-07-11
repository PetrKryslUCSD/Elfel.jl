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
FESpaces.ndofsperel
FESpaces.numberfreedofs!
FESpaces.numberdatadofs!
FESpaces.ndofs
FESpaces.nunknowns
FESpaces.highestfreedofnum
FESpaces.highestdatadofnum
FESpaces.numberdofs!
FESpaces.setebc!
FESpaces.gathersysvec!
FESpaces.gathersysvec!
FESpaces.scattersysvec!
FESpaces.scattersysvec!
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
FEIterators.ndofsperel
FEIterators.eldofs
FEIterators.elnodes
FEIterators.eldofentmdims
FEIterators.eldofcomps
FEIterators.jacjac(it::FEIterators.FEIterator, qpit::QPIterators.QPIterator)
```

## Quadrature-point iterators

```@docs
QPIterators.QPIterator
Base.iterate
QPIterators.bfun
QPIterators.bfungradpar
QPIterators.bfungrad
QPIterators.weight
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

