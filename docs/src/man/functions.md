# Functions

## Reference shapes

```@meta
CurrentModule = Elfel.RefShapes
```

```@docs
Elfel.RefShapes.manifdim
Elfel.RefShapes.manifdimv
Elfel.RefShapes.quadrature
```


## Elements

```@meta
CurrentModule = Elfel.FElements
```

```@docs
Elfel.FElements.shapedesc
Elfel.FElements.refshape
Elfel.FElements.nfeatofdim
Elfel.FElements.ndofperfeat
Elfel.FElements.ndofsperel
Elfel.FElements.manifdim
Elfel.FElements.Jacobian
Elfel.FElements.jacjac
Elfel.FElements.bfun
Elfel.FElements.bfungradpar
Elfel.FElements.FEH1_L2
Elfel.FElements.FEH1_T3
Elfel.FElements.FEH1_Q4
Elfel.FElements.FEH1_T6
Elfel.FElements.FEH1_T3_BUBBLE
```

## Fields

```@docs
Elfel.FEFields.doftype
Elfel.FEFields.dofnumtype
Elfel.FEFields.nterms
Elfel.FEFields.ndofsperterm
Elfel.FEFields.ndofs
Elfel.FEFields.setebc!
Elfel.FEFields.numberfreedofs!
Elfel.FEFields.numberdatadofs!
Elfel.FEFields.freedofnums
Elfel.FEFields.datadofnums
Elfel.FEFields.highestfreedofnum
Elfel.FEFields.highestdatadofnum
Elfel.FEFields.gathersysvec!
Elfel.FEFields.scattersysvec!
```

## Spaces

```@meta
CurrentModule = Elfel.FESpaces
```


```@docs
Elfel.FESpaces.doftype
Elfel.FESpaces.edofmdim
Elfel.FESpaces.edofbfnum
Elfel.FESpaces.edofcompnt
Elfel.FESpaces.ndofsperel
Elfel.FESpaces.numberfreedofs!
Elfel.FESpaces.numberdatadofs!
Elfel.FESpaces.ndofs
Elfel.FESpaces.nunknowns
Elfel.FESpaces.highestfreedofnum
Elfel.FESpaces.highestdatadofnum
Elfel.FESpaces.numberdofs!
Elfel.FESpaces.setebc!
Elfel.FESpaces.gathersysvec!(v, fesp::FESpace)
Elfel.FESpaces.gathersysvec!(v, fesp::AbstractVector) 
Elfel.FESpaces.scattersysvec!(fesp::FESpace, v)
Elfel.FESpaces.scattersysvec!(fesp::AbstractVector, v) 
Elfel.FESpaces.makeattribute
```

## Finite element iterators

```@docs
Base.iterate
Elfel.FEIterators.ndofsperel
Elfel.FEIterators.eldofs
Elfel.FEIterators.elnodes
Elfel.FEIterators.eldofentmdims
Elfel.FEIterators.eldofcomps
Elfel.FEIterators.jacjac
```

## Quadrature-point iterators

```@docs
Base.iterate
Elfel.QPIterators.bfun
Elfel.QPIterators.bfungradpar
Elfel.QPIterators.bfungrad
Elfel.QPIterators.weight
```

## Assemblers

```@docs
Elfel.Assemblers.SysmatAssemblerSparse
Elfel.Assemblers.start!
Elfel.Assemblers.assemble!
Elfel.Assemblers.finish!
Elfel.Assemblers.SysvecAssembler
```

## Local Assemblers

```@docs
Elfel.LocalAssemblers.LocalMatrixAssembler
Elfel.LocalAssemblers.LocalVectorAssembler
Base.size
Base.getindex
Base.setindex!
Elfel.LocalAssemblers.init!
```

## Index

```@index
```

