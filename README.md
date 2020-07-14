[![Project Status: Active – The project is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://img.shields.io/travis/PetrKryslUCSD/Elfel.jl/master.svg?label=Linux+MacOSX+Windows)](https://travis-ci.org/PetrKryslUCSD/Elfel.jl)
[![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/Elfel.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/Elfel.jl?branch=master)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/Elfel.jl/dev)

# Elfel.jl: Extensible library for finite element methods

**elfel**<br>
G. *adjective*. different, like something else, unique, unparalled

## Introduction

This package provides support for the development of Finite Element Method applications, especially in continuum mechanics. Mixed methods with cooperating finite element spaces are supported. High performance is one of the focus points. Check out the [tutorials](https://petrkryslucsd.github.io/Elfel.jl/dev/#Tutorials) to get a feel for the functionality of this package.

## Dependencies

The library relies on mesh support from the following suite of small mesh-management packages: [`MeshCore`](https://github.com/PetrKryslUCSD/MeshCore.jl), [`MeshSteward`](https://github.com/PetrKryslUCSD/MeshSteward.jl).

## News

- 07/11/2020: Tutorial structure added.
- 07/06/2020: Exports have been added to facilitate use of the library.
- 07/05/2020: Vector finite element spaces tested.
- 06/26/2020: L2 elements  implemented, Stokes problem example.
- 06/19/2020: Example of a mixed method for the discrete Stokes problem added.
- 05/25/2020: Firmed up the concept of iterators for access to element and quadrature point data.


## Usage

Clone the repository, and execute
```
using Pkg; Pkg.activate("."); Pkg.instantiate()
```
in the `Elfel` folder.

The user can either use/import individual functions from `Elfel` like so:
```
using Elfel.FElements: FEH1_T6, FEH1_T3, refshape, Jacobian
```
or all exported symbols maybe made available in the user's context as
```
using Elfel.Exports
```

The `examples` folder can be explored by simply running the files with `include()`.

