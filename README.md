# Elfel.jl: Extensible library for finite element methods

**elfel**<br>
G. *adjective*. different, like something else, unique, unparalled

## Introduction

This package provides support for the development of Finite Element Method applications, especially in continuum mechanics.

## Dependencies

The library relies on mesh support from the following suite of small mesh-management packages:
```
github.com/PetrKryslUCSD/MeshCore.jl
github.com/PetrKryslUCSD/MeshSteward.jl     
```

## News

- 05/25/2020: Firmed up the concept of iterators for access to element and quadrature point data.


## Usage

Clone the repository, and execute
```
using Pkg; Pkg.activate("."); Pkg.instantiate()
```
in the `Elfel` folder.