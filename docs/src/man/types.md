# Types

## Reference shapes

```@meta
CurrentModule = Elfel.RefShapes
```


```@docs
AbstractRefShape
RefShapePoint
RefShapeInterval
RefShapeSquare
RefShapeCube
RefShapeTriangle
RefShapeTetrahedron
IntegRule
```


## Elements

```@meta
CurrentModule = Elfel.FElements
```

```@docs
FE{RS, SD}
FEData{SD}
```

## Fields

```@meta
CurrentModule = Elfel.FEFields
```

```@docs
FEField{N, T, IT}
```

## Spaces

```@meta
CurrentModule = Elfel.FESpaces
```


```@docs
FESpace{FET, T}
```

## Finite element iterators

```@meta
CurrentModule = Elfel.FEIterators
```

```@docs
FEIterator
```

## Quadrature-point iterators

```@meta
CurrentModule = Elfel.QPIterators
```

```@docs
QPIterator
```

## Assemblers

```@meta
CurrentModule = Elfel.Assemblers
```

```@docs
AbstractSysmatAssembler
SysmatAssemblerSparse
SysvecAssembler
```

## Local Assemblers

```@meta
CurrentModule = Elfel.LocalAssemblers
```

```@docs
LocalMatrixAssembler
LocalVectorAssembler
```

## Index

```@index
```

