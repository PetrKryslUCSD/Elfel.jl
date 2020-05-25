module FElements

using StaticArrays
using LinearAlgebra
using MeshCore
using MeshCore: manifdim
using ..RefShapes: RefShapePoint, RefShapeInterval, RefShapeTriangle, RefShapeTetrahedron, RefShapeSquare, RefShapeCube


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
    ndofsperelem(fe::FE{RS, SD}) where {RS, SD}

Provide the number of nodes per element.
"""
function ndofsperelem(fe::FE{RS, SD}) where {RS, SD}
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
    gradN!(::Val{1}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates, one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(::Val{1}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}
    r = 1.0 / redJ[1, 1]
    for r in 1:size(gradN, 1)
        gradN[r, 1] =  gradNparams[r, 1] * r
    end
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
    gradN!(::Val{2}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates,
  one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(::Val{2}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}
    # This is the unrolled version that avoids allocation of a 2 x 2 matrix
    invdet = 1.0/(redJ[1, 1]*redJ[2, 2] - redJ[1, 2]*redJ[2, 1]);
    invredJ11 =  (redJ[2, 2])*invdet;
    invredJ12 = -(redJ[1, 2])*invdet;
    invredJ21 = -(redJ[2, 1])*invdet;
    invredJ22 =  (redJ[1, 1])*invdet;
    @assert size(gradN, 1)==size(gradNparams, 1)
    @inbounds for r in 1:size(gradN, 1)
        gradN[r, 1] = gradNparams[r, 1]*invredJ11 + gradNparams[r, 2]*invredJ21;
        gradN[r, 2] = gradNparams[r, 1]*invredJ12 + gradNparams[r, 2]*invredJ22;
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

"""
    gradN!(::Val{3}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}

Compute the gradient of the basis functions with the respect to
the "reduced" spatial coordinates.

- `gradN`= output,  matrix of gradients,  one per row
- `gradNparams`= matrix of gradients with respect to parametric coordinates,
  one per row
- `redJ`= reduced Jacobian matrix `redJ=transpose(Rm)*J`
"""
function gradN!(::Val{3}, gradN::T1, gradNparams::T2, redJ::T3) where {T1, T2, T3}
    invdet = 1.0 / ( +redJ[1, 1]*(redJ[2, 2]*redJ[3, 3]-redJ[3, 2]*redJ[2, 3])
                    -redJ[1, 2]*(redJ[2, 1]*redJ[3, 3]-redJ[2, 3]*redJ[3, 1])
                    +redJ[1, 3]*(redJ[2, 1]*redJ[3, 2]-redJ[2, 2]*redJ[3, 1]) );
    # This is the unrolled version that avoids allocation of a 3 x 3 matrix
    invredJ11 =  (redJ[2, 2]*redJ[3, 3]-redJ[3, 2]*redJ[2, 3])*invdet;
    invredJ12 = -(redJ[1, 2]*redJ[3, 3]-redJ[1, 3]*redJ[3, 2])*invdet;
    invredJ13 =  (redJ[1, 2]*redJ[2, 3]-redJ[1, 3]*redJ[2, 2])*invdet;
    invredJ21 = -(redJ[2, 1]*redJ[3, 3]-redJ[2, 3]*redJ[3, 1])*invdet;
    invredJ22 =  (redJ[1, 1]*redJ[3, 3]-redJ[1, 3]*redJ[3, 1])*invdet;
    invredJ23 = -(redJ[1, 1]*redJ[2, 3]-redJ[2, 1]*redJ[1, 3])*invdet;
    invredJ31 =  (redJ[2, 1]*redJ[3, 2]-redJ[3, 1]*redJ[2, 2])*invdet;
    invredJ32 = -(redJ[1, 1]*redJ[3, 2]-redJ[3, 1]*redJ[1, 2])*invdet;
    invredJ33 =  (redJ[1, 1]*redJ[2, 2]-redJ[2, 1]*redJ[1, 2])*invdet;
    @assert size(gradN, 1)==size(gradNparams, 1)
    @inbounds for r in 1:size(gradN, 1)
        gradN[r, 1] = gradNparams[r, 1]*invredJ11 + gradNparams[r, 2]*invredJ21 + gradNparams[r, 3]*invredJ31;
        gradN[r, 2] = gradNparams[r, 1]*invredJ12 + gradNparams[r, 2]*invredJ22 + gradNparams[r, 3]*invredJ32;
        gradN[r, 3] = gradNparams[r, 1]*invredJ13 + gradNparams[r, 2]*invredJ23 + gradNparams[r, 3]*invredJ33;
    end
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

# T3 ==================================================================
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
FEH1_Q4_TYPE = FE{RefShapeTriangle, typeof(MeshCore.Q4)}
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
