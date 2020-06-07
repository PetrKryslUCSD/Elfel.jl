module FElements

using StaticArrays
using LinearAlgebra
using MeshCore
using MeshCore: manifdim
using ..RefShapes: RefShapePoint, RefShapeInterval, RefShapeTriangle, RefShapeTetrahedron, RefShapeSquare, RefShapeCube, manifdimv


"""
    FE{RS, SD} 

Type of a finite element set. 

- `RS` = reference shape,
- `SD` = shape descriptor.
"""
struct FE{RS, SD} 
    sd::SD
    ndof::SVector{4, Int64}
end

"""
    shapedesc(fe::FE{RS, SD}) where {RS, SD}

Topological shape description.

Refer to the MeshCore library.
"""
shapedesc(fe::FE{RS, SD}) where {RS, SD} = fe.sd

"""
    refshape(fe::FE{RS, SD}) where {RS, SD}

Reference shape.
"""
refshape(fe::FE{RS, SD}) where {RS, SD} = RS

"""
    nfeatofdim(fe::FE{RS, SD}, m) where {RS, SD}

Number of features of manifold dimension `m`.
"""
nfeatofdim(fe::FE{RS, SD}, m) where {RS, SD} = MeshCore.nfeatofdim(fe.sd, m)

"""
    ndofsperfeat(fe::FE{RS, SD}, m) where {RS, SD}

Number of degrees of freedom per feature of manifold dimension `m`. 
"""
ndofsperfeat(fe::FE{RS, SD}, m) where {RS, SD} = fe.ndof[m+1]

"""
    ndofsperel(fe::FE{RS, SD}) where {RS, SD}

Provide the number of degrees of freedom per element.
Enumerate all features of all manifold dimensions, and for each feature multiply
by the number of degrees of freedom per feature.
"""
function ndofsperel(fe::FE{RS, SD}) where {RS, SD}
    md = manifdim(fe.sd)
    n = 0
    for m in 0:1:md
        n = n + nfeatofdim(fe, m) * ndofsperfeat(fe, m)
    end
    return n
end

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
    NBFPE = length(gradNpar)
    j = 1
    J = locs[conn[j]] * gradNpar[j]
    @inbounds for j in 2:NBFPE
        J += locs[conn[j]] * gradNpar[j]
    end
    return J
end

"""
    jacjac(fe::FE{RS, SD}, locs, gradNpar) where {RS, SD}

Compute the Jacobian matrix and the Jacobian determinant.

This is the generic version suitable for isoparametric elements.
"""
function jacjac(fe::FE{RS, SD}, locs, nodes, gradNpar) where {RS, SD}
    Jac = _jac(locs, nodes, gradNpar)
    return (Jac, Jacobian(manifdimv(refshape(fe)), Jac))
end

# L2 ==================================================================
FEH1_L2_TYPE = FE{RefShapeInterval, typeof(MeshCore.L2)}
FEH1_L2(NDN) = FEH1_L2_TYPE(MeshCore.L2, SVector{4}([NDN, 0, 0, 0]))

function bfun(self::FEH1_L2_TYPE,  param_coords) 
    return SVector{2}([(1. - param_coords[1]); (1. + param_coords[1])] / 2.0)
end

function bfungradpar(self::FEH1_L2_TYPE,  param_coords) 
    g = reshape([-1.0; +1.0]/2.0, 2, 1)
    return [SVector{1}(g[idx, :])' for idx in 1:size(g, 1)]
end

# T3 
FEH1_T3_TYPE = FE{RefShapeTriangle, typeof(MeshCore.T3)}
FEH1_T3(NDN) = FEH1_T3_TYPE(MeshCore.T3, SVector{4}([NDN, 0, 0, 0]))

function bfun(self::FEH1_T3_TYPE,  param_coords) 
    return SVector{3}([(1 - param_coords[1] - param_coords[2]); param_coords[1]; param_coords[2]])
end

function bfungradpar(self::FEH1_T3_TYPE,  param_coords)
    g = [-1. -1.;  +1.  0.;  0. +1.]
    return [SVector{2}(g[idx, :])' for idx in 1:size(g, 1)]
end

# Q4 ==================================================================
FEH1_Q4_TYPE = FE{RefShapeSquare, typeof(MeshCore.Q4)}
FEH1_Q4(NDN) = FEH1_Q4_TYPE(MeshCore.Q4, SVector{4}([NDN, 0, 0, 0]))

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

end # module
