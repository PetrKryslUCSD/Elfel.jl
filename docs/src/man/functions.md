# Functions

## Reference shapes

```@meta
CurrentModule = Elfel.RefShapes
```


```@docs
manifdim
manifdimv
quadrature
```


## Elements

```@meta
CurrentModule = Elfel.FElements
```

```@docs
shapedesc
refshape
nfeatofdim
ndofperfeat
ndofsperel(fe::FE{RS, SD}) where {RS, SD}
manifdim
Jacobian
jacjac
bfun
bfungradpar
FEH1_L2
FEH1_T3
FEH1_Q4
FEH1_T6
FEH1_T3_BUBBLE
```

## Fields

```@meta
CurrentModule = Elfel.FEFields
```

```@docs
doftype
dofnumtype
nterms
ndofsperterm
ndofs
setebc!
numberfreedofs!
numberdatadofs!
freedofnums
datadofnums
highestfreedofnum
highestdatadofnum
gathersysvec!
scattersysvec!
```

## Spaces

```@meta
CurrentModule = Elfel.FESpaces
```


```@docs
doftype(fesp::FESpace{FET, T}) where {FET, T}
edofmdim
edofbfnum
edofcompnt
ndofsperel(fesp::FES)  where {FES<:FESpace}
numberfreedofs!(fesp::FES, firstnum = 1)  where {FES<:FESpace}
numberdatadofs!(fesp::FES, firstnum = 0)  where {FES<:FESpace}
ndofs(fesp::FES)  where {FES<:FESpace}
nunknowns(fesp::FES)  where {FES<:FESpace}
highestfreedofnum(fesp::FES)  where {FES<:FESpace}
highestdatadofnum(fesp::FES)  where {FES<:FESpace}
numberdofs!
setebc!(fesp::FESpace, mid, eid, comp, val::T) where {T}
gathersysvec!(v, fesp::FESpace)
gathersysvec!(v, fesp::AbstractVector) 
scattersysvec!(fesp::FESpace, v)
scattersysvec!(fesp::AbstractVector, v) 
makeattribute
```

## Finite element iterators

```@meta
CurrentModule = Elfel.FEIterators
```

```@docs
Base.iterate
ndofsperel
eldofs
elnodes
eldofentmdims
eldofcomps
jacjac(it::FEIterator, qpit::QPIterator)
```

## Quadrature-point iterators

```@meta
CurrentModule = Elfel.QPIterators
```

```@docs
Base.iterate
bfun
bfungradpar
bfungrad
weight
```

## Assemblers

```@meta
CurrentModule = Elfel.Assemblers
```

```@docs
SysmatAssemblerSparse
start!
assemble!
finish!
SysvecAssembler
```

## Local Assemblers

```@meta
CurrentModule = Elfel.LocalAssemblers
```

```@docs
LocalMatrixAssembler
LocalVectorAssembler
Base.size
Base.getindex
Base.setindex!
init!
```

## Index

```@index
```

