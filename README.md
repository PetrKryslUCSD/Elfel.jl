[![Project Status: Active – The project is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://img.shields.io/travis/PetrKryslUCSD/Elfel.jl/master.svg?label=Linux+MacOSX+Windows)](https://travis-ci.org/PetrKryslUCSD/Elfel.jl)
[![Coverage Status](https://coveralls.io/repos/github/PetrKryslUCSD/Elfel.jl/badge.svg?branch=master)](https://coveralls.io/github/PetrKryslUCSD/Elfel.jl?branch=master)
[![Latest documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://petrkryslucsd.github.io/Elfel.jl/dev)

# Elfel.jl: Extensible library for finite element methods

**elfel**<br>
G. *adjective*. different, like something else, unique, unparalled

## Introduction
## 
This package provides support for the development of Finite Element Method applications, especially in continuum mechanics.

## Dependencies

The library relies on mesh support from the following suite of small mesh-management packages:
```
github.com/PetrKryslUCSD/MeshCore.jl
github.com/PetrKryslUCSD/MeshSteward.jl     
```

## News

- 06/19/2020: Example of a mixed method for the discrete Stokes problem added.
- 05/25/2020: Firmed up the concept of iterators for access to element and quadrature point data.


## Usage

Clone the repository, and execute
```
using Pkg; Pkg.activate("."); Pkg.instantiate()
```
in the `Elfel` folder.