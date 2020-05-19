module mfexp1
using Elfel
using Elfel.RefShapes: RefShapeTriangle, manifdim, RefShapeInterval
using Elfel.FElements: FE, nodesperelem, refshape
using Elfel.FElements: bfun, bfundpar, nbasisfuns
using Elfel.Fields: FEField, ndofsperentity, nentities, numberdofs!, setebc!, gathersysvec
using Elfel.FESpaces: FESpace
using Elfel.FEExpansions: FEExpansion
using MeshKeeper: Mesh, load
using Test
function test()
    e = FE{RefShapeTriangle, 3, 1}()

    fesp = FESpace((e, 1))
    mesh = load(Mesh(), "mt3gen3.mesh")
    fex = FEExpansion(mesh, fesp)
end
end
using .mfexp1
mfexp1.test()
