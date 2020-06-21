module FElements

using StaticArrays
using LinearAlgebra
using MeshCore
import MeshCore: manifdim
using ..RefShapes: RefShapePoint, RefShapeInterval, RefShapeTriangle, RefShapeTetrahedron, RefShapeSquare, RefShapeCube, manifdimv

"""
    FE{RS, SD}

Abstract type of finite element, parameterized by
- `RS`: type of reference shape of the element (triangle, square, ...), and
- `SD`: shape descriptor; refer to the package `MeshCore`.
"""
abstract type FE{RS, SD} end

"""
    FEData{SD} 

Type of a finite element data. 

Parameterized by
- `SD` = shape descriptor; refer to the package `MeshCore`.
"""
struct FEData{SD} 
    sd::SD
    _ndofperfeat::SVector{4, Int64}
end

"""
    shapedesc(fe::FE{RS, SD}) where {RS, SD}

Topological shape description.

Refer to the MeshCore library.
"""
shapedesc(fe::FE{RS, SD}) where {RS, SD} = fe.data.sd

"""
    refshape(fe::FE{RS, SD}) where {RS, SD}

Reference shape.
"""
refshape(fe::FE{RS, SD}) where {RS, SD} = RS

"""
    nfeatofdim(fe::FE{RS, SD}, m) where {RS, SD}

Number of features of manifold dimension `m`. Note that `0 <= m <= 3`.
"""
nfeatofdim(fe::FE{RS, SD}, m) where {RS, SD} = MeshCore.nfeatofdim(fe.data.sd, m)

"""
    feathasdof(fe::FE{RS, SD}, m) where {RS, SD}

How many degrees of freedom are attached to a the feature of manifold dimension
`m`?

Note that `0 <= m <= 3`.
"""
ndofperfeat(fe::FE{RS, SD}, m) where {RS, SD} = fe.data._ndofperfeat[m+1]

"""
    ndofsperel(fe::FE{RS, SD}) where {RS, SD}

Provide the number of degrees of freedom per element.

Enumerate all features of all manifold dimensions, and for each feature multiply
by the number of degrees of freedom per feature. The assumption is that this is
a *scalar* finite element.
"""
function ndofsperel(fe::FE{RS, SD}) where {RS, SD}
    md = manifdim(fe.data.sd)
    n = 0
    for m in 0:1:md
        n = n + nfeatofdim(fe, m) * ndofperfeat(fe, m)
    end
    return n
end

"""
    manifdim(fe::FE{RS, SD}) where {RS, SD}

Get the manifold dimension of the finite element.
"""
manifdim(fe::FE{RS, SD}) where {RS, SD} = MeshCore.manifdim(fe.data.sd)

"""
    Jacobian(::Val{0}, J::T) where {T}

Evaluate the point Jacobian.

- `J` = Jacobian matrix, which isn't really defined well for a 0-manifold.
"""
function Jacobian(::Val{0}, J::T) where {T}
    return 1.0;
end

"""
    Jacobian(::Val{1}, J::T) where {T}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(::Val{1}, J::T) where {T}
    sdim,  ntan = size(J);
    @assert ntan == 1 "Expected number of tangent vectors: 1"
    return norm(J)
end

"""
    Jacobian(::Val{2}, J::T) where {T}

Evaluate the curve Jacobian.

- `J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(::Val{2}, J::T) where {T}
    sdim,  ntan = size(J);
    @assert ntan == 2 "Expected number of tangent vectors: 2"
    if sdim == ntan
        @inbounds Jac = (J[1, 1]*J[2, 2] - J[2, 1]*J[1, 2])
        return Jac;# is det(J);% Compute the Jacobian
    else
        return norm(cross(J[:, 1], J[:, 2]));
    end
end

"""
    Jacobian(self::T, J::FFltMat)::FFlt where {T<:AbstractFESet3Manifold}

Evaluate the volume Jacobian.

`J` = Jacobian matrix, columns are tangent to parametric coordinates curves.
"""
function Jacobian(::Val{3}, J::T) where {T}
    sdim,  ntan = size(J);
    @assert (ntan == 3) && (sdim == 3) "Expected number of tangent vectors: 3"
    #Jac = det(J);# Compute the Jacobian
    # The unrolled version
    return (+J[1, 1]*(J[2, 2]*J[3, 3]-J[3, 2]*J[2, 3])
    -J[1, 2]*(J[2, 1]*J[3, 3]-J[2, 3]*J[3, 1])
    +J[1, 3]*(J[2, 1]*J[3, 2]-J[2, 2]*J[3, 1]) )
end

function _jac(locs, conn, gradNpar)
    NBFPE = length(conn)
    j = 1
    J = locs[conn[j]] * gradNpar[j]
    @inbounds for j in 2:NBFPE
        J += locs[conn[j]] * gradNpar[j]
    end
    return J
end

"""
    jacjac(fe::FE{RS, SD}, locs, nodes, gradNpar) where {RS, SD}

Compute the Jacobian matrix and the Jacobian determinant.

This is the generic version suitable for isoparametric elements.
"""
function jacjac(fe::FE{RS, SD}, locs, nodes, gradNpar) where {RS, SD}
    Jac = _jac(locs, nodes, gradNpar)
    return (Jac, Jacobian(manifdimv(refshape(fe)), Jac))
end

"""
    bfun(self::FESUBT,  param_coords)  where {FESUBT<:FE{RS, SD}}

Evaluate the basis functions for all degrees of freedom of the scalar finite
element at the parametric coordinates. Return a vector of the values.
"""
function bfun(self::FESUBT,  param_coords)  where {FESUBT<:FE{RS, SD}}
end

"""
    bfungradpar(self::FESUBT,  param_coords)  where {FESUBT<:FE{RS, SD}}

Evaluate the gradients of the basis functions for all degrees of freedom of
the scalar finite element with respect to the parametric coordinates, at the
parametric coordinates given. Return a vector of the gradients.
"""
function bfungradpar(self::FESUBT,  param_coords)  where {FESUBT<:FE{RS, SD}}
end

# L2 ==================================================================
# Linear two-node element. Only nodal basis functions.
struct FEH1_L2_Type{RS, SD} <: FE{RS, SD}
    data::FEData{SD}
end
FEH1_L2_TYPE = FEH1_L2_Type{RefShapeInterval, typeof(MeshCore.L2)}

"""
    FEH1_L2()

Construct an H1 finite element of the type L2.

L2 is two-node linear segment element.
"""
FEH1_L2() = FEH1_L2_TYPE(FEData(MeshCore.L2, SVector{4}([1, 0, 0, 0])))

function bfun(self::FEH1_L2_TYPE,  param_coords) 
    return SVector{2}([(1. - param_coords[1]); (1. + param_coords[1])] / 2.0)
end

function bfungradpar(self::FEH1_L2_TYPE,  param_coords) 
    g = reshape([-1.0; +1.0]/2.0, 2, 1)
    return [SVector{1}(g[idx, :])' for idx in 1:size(g, 1)]
end

# T3 ==================================================================
# Linear triangular element. Only nodal basis functions.
struct FEH1_T3_Type{RS, SD} <: FE{RS, SD}
    data::FEData{SD}
end
FEH1_T3_TYPE = FEH1_T3_Type{RefShapeTriangle, typeof(MeshCore.T3)}

"""
    FEH1_T3()

Construct an H1 finite element of the type T3.

T3 is 3-node linear triangle element.
"""
FEH1_T3() = FEH1_T3_TYPE(FEData(MeshCore.T3, SVector{4}([1, 0, 0, 0])))

function bfun(self::FEH1_T3_TYPE,  param_coords) 
    return SVector{3}([(1 - param_coords[1] - param_coords[2]); param_coords[1]; param_coords[2]])
end

function bfungradpar(self::FEH1_T3_TYPE,  param_coords)
    g = [-1. -1.;  +1.  0.;  0. +1.]
    return [SVector{2}(g[idx, :])' for idx in 1:size(g, 1)]
end

# T6 ==================================================================
# Quadratic triangular element. Only nodal basis functions.
struct FEH1_T6_Type{RS, SD} <: FE{RS, SD}
    data::FEData{SD}
end
FEH1_T6_TYPE = FEH1_T6_Type{RefShapeTriangle, typeof(MeshCore.T6)}

"""
    FEH1_T6()

Construct an H1 finite element of the type T6.

T6 is 6-node quadratic triangle element.
"""
FEH1_T6() = FEH1_T6_TYPE(FEData(MeshCore.T6, SVector{4}([1, 0, 0, 0])))

function bfun(self::FEH1_T6_TYPE,  param_coords) 
    r=param_coords[1];
    s=param_coords[2];
    t = 1. - r - s;
    val = [t * (t + t - 1);
           r * (r + r - 1);
           s * (s + s - 1);
           4 * r * t;
           4 * r * s;
           4 * s * t];
    return SVector{6}(val)
end

function bfungradpar(self::FEH1_T6_TYPE,  param_coords)
    r =param_coords[1];
    s =param_coords[2];
    t = 1. - r - s;
    val = [-3+4*r+4*s  -3+4*r+4*s;
           4*r-1  0.0;
           0.0 4*s-1;
           4-8*r-4*s  -4*r;
           4*s  4*r;
           -4*s  4-4*r-8*s];
    return [SVector{2}(val[idx, :])' for idx in 1:size(val, 1)]
end

# Q4 ==================================================================
# Linear quadrilateral element. Only nodal basis functions.
struct FEH1_Q4_Type{RS, SD} <: FE{RS, SD}
    data::FEData{SD}
end
FEH1_Q4_TYPE = FEH1_Q4_Type{RefShapeSquare, typeof(MeshCore.Q4)}

"""
    FEH1_Q4()

Construct an H1 finite element of the type Q4.

Q4 is 4-node linear quadrilateral element.
"""
FEH1_Q4() = FEH1_Q4_TYPE(FEData(MeshCore.Q4, SVector{4}([1, 0, 0, 0])))

function bfun(self::FEH1_Q4_TYPE,  param_coords) 
	val = [0.25 * (1. - param_coords[1]) * (1. - param_coords[2]);
	       0.25 * (1. + param_coords[1]) * (1. - param_coords[2]);
	       0.25 * (1. + param_coords[1]) * (1. + param_coords[2]);
	       0.25 * (1. - param_coords[1]) * (1. + param_coords[2])];
    return SVector{4}(val)
end

function bfungradpar(self::FEH1_Q4_TYPE,  param_coords) 
    g =   [-(1. - param_coords[2])*0.25 -(1. - param_coords[1])*0.25;
            (1. - param_coords[2])*0.25 -(1. + param_coords[1])*0.25;
            (1. + param_coords[2])*0.25 (1. + param_coords[1])*0.25;
           -(1. + param_coords[2])*0.25 (1. - param_coords[1])*0.25];
    return [SVector{2}(g[idx, :])' for idx in 1:size(g, 1)]
end

# T3-BUBBLE ==================================================================
# Linear triangular element with a cubic interior bubble. 
struct FEH1_T3_BUBBLE_Type{RS, SD} <: FE{RS, SD}
    data::FEData{SD}
end
FEH1_T3_BUBBLE_TYPE = FEH1_T3_BUBBLE_Type{RefShapeTriangle, typeof(MeshCore.T3)}

"""
    FEH1_T3_BUBBLE()

Construct an H1 finite element of the type T3 with a cubic bubble.

T3 is 3-node linear triangle element with a cubic bubble. It has the usual
nodal basis functions associated with the vertices, and cubic bubble
associated with the element itself.
"""
FEH1_T3_BUBBLE() = FEH1_T3_BUBBLE_TYPE(FEData(MeshCore.T3, SVector{4}([1, 0, 1, 0])))

function bfun(self::FEH1_T3_BUBBLE_TYPE,  param_coords) 
    xi, eta = param_coords
    return SVector{4}([(1 - xi - eta); 
        xi; 
        eta; 
        (1 - xi - eta) * xi * eta])
end

function bfungradpar(self::FEH1_T3_BUBBLE_TYPE,  param_coords)
    xi, eta = param_coords
    g = [-1. -1.;  
    +1.  0.;  
    0. +1.;
    (-xi * eta + (1 - xi - eta) * eta) (-xi * eta + (1 - xi - eta) * xi)]
    return [SVector{2}(g[idx, :])' for idx in 1:size(g, 1)]
end

end # module
