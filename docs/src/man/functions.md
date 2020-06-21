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

## Fields

```@docs
FEFields.doftype
FEFields.dofnumtype
FEFields.nterms
FEFields.ndofsperterm
FEFields.ndofs
FEFields.setebc!
FEFields.numberfreedofs!
FEFields.numberdatadofs!
FEFields.freedofnums
FEFields.datadofnums
FEFields.highestfreedofnum
FEFields.highestdatadofnum
FEFields.gathersysvec!
FEFields.scattersysvec!
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
FESpaces.gathersysvec!(v, fesp)
FESpaces.gathersysvec!(v, fesp::AbstractVector) 
FESpaces.scattersysvec!(fesp, v)
FESpaces.scattersysvec!(fesp::AbstractVector, v) 
FESpaces.makeattribute
```

## Finite element iterators

```@docs
Base.iterate
FEIterators.ndofsperel
FEIterators.eldofs
FEIterators.elnodes
FEIterators.eldofentmdims
FEIterators.eldofcomps
FEIterators.jacjac
```

## Quadrature-point iterators

```@docs
Base.iterate
QPIterators.bfun(it::QPIterators.QPIterator)
QPIterators.bfungradpar(it::QPIterators.QPIterator)
QPIterators.bfungrad
QPIterators.weight
```

## Assemblers

```@docs
Assemblers.SysmatAssemblerSparse
Assemblers.start!
Assemblers.assemble!
Assemblers.finish!
Assemblers.SysvecAssembler
```

## Local Assemblers

```@docs
LocalAssemblers.LocalMatrixAssembler
LocalAssemblers.LocalVectorAssembler
Base.size
Base.getindex
Base.setindex!
LocalAssemblers.init!
```

## Index

```@index
```

