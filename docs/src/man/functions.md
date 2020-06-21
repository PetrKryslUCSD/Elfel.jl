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
FESpaces.ndofsperel(fesp::FESpaces.FESpace)
FESpaces.numberfreedofs!(fesp::FESpaces.FESpace, firstnum = 1)
FESpaces.numberdatadofs!(fesp::FESpaces.FESpace, firstnum = 0)
FESpaces.ndofs
FESpaces.nunknowns
FESpaces.highestfreedofnum(fesp::FESpaces.FESpace)
FESpaces.highestdatadofnum(fesp::FESpaces.FESpace)
FESpaces.numberdofs!
FESpaces.setebc!(fesp::FESpaces.FESpace, mid, eid, comp, val)
FESpaces.gathersysvec!(v, fesp::FESpaces.FESpace)
FESpaces.gathersysvec!(v, fesp::AbstractVector) 
scattersysvec!(fesp::FESpaces.FESpace, v)
FESpaces.scattersysvec!(fesp::AbstractVector, v) 
FESpaces.makeattribute
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

