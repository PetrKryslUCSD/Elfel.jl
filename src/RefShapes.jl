module RefShapes

using LinearAlgebra

"""
    AbstractRefShape{MANIFDIM}

Abstract type of a reference shape.
"""
abstract type AbstractRefShape{MANIFDIM} end

"""
    RefShapePoint <: AbstractRefShape{0}

Type of a reference shape for a zero-dimensional manifold (point).
"""
struct RefShapePoint <: AbstractRefShape{0} end

"""
    RefShapeInterval <: AbstractRefShape{1}

Type of a reference shape for a 1-dimensional manifold (curve).
"""
struct RefShapeInterval <: AbstractRefShape{1} end

"""
    RefShapeSquare <: AbstractRefShape{2}

Type of a logically rectangular reference shape for a 2-dimensional manifold 
(surface).
"""
struct RefShapeSquare <: AbstractRefShape{2} end

"""
    RefShapeCube <: AbstractRefShape{3}

Type of a reference shape for a 3-dimensional manifold (solid) bounded by six quadrilaterals.
"""
struct RefShapeCube <: AbstractRefShape{3} end

"""
    RefShapeTriangle <: AbstractRefShape{2}

Type of a logically triangular reference shape for a 2-dimensional manifold 
(surface).
"""
struct RefShapeTriangle <: AbstractRefShape{2} end

"""
    RefShapeTetrahedron <: AbstractRefShape{3}

Type of a reference shape for a 3-dimensional manifold (solid) bounded by 4 triangles.
"""
struct RefShapeTetrahedron <: AbstractRefShape{3} end

"""
    manifdim(rs)

Get the manifold dimension of the reference shape.
"""
manifdim(::Type{T}) where {T<:AbstractRefShape{MANIFDIM}} where {MANIFDIM} = MANIFDIM

"""
    manifdimv(::Type{T}) where {T<:AbstractRefShape{MANIFDIM}} where {MANIFDIM}

Get the manifold dimension of the reference shape as `Val`. 
"""
manifdimv(::Type{T}) where {T<:AbstractRefShape{MANIFDIM}} where {MANIFDIM} = Val(MANIFDIM)

"""
    IntegRule

Type for integration rule.
"""
struct IntegRule 
    npts::Int64
    param_coords::Array{Float64, 2}
    weights::Array{Float64, 1}
end

npts(qr) = qr.npts
param_coords(qr) = qr.param_coords
weights(qr) = qr.weights

function _gauss1(order)
    function opengauss(N, a = -1.0, b = +1.0) 
        F = eigen(SymTridiagonal(zeros(N), [n/sqrt(4n^2 - 1) for n = 1:N-1])) 
        return [(F.values[i]+1)*(b-a)/2 + a for i = 1:N], [2*F.vectors[1, i]^2 for i = 1:N]*(b-a)/2 
    end
    
    if     (order==1)
        param_coords = vec([ 0.0 ]);
        weights  = vec([ 2.0 ]);
    elseif (order==2)
        param_coords = vec([ -0.577350269189626 0.577350269189626 ]);
        weights  = vec([ 1.0 1.0 ]);
    elseif (order==3)
        param_coords = vec([ -0.774596669241483  0.0  0.774596669241483 ]);
        weights  = vec([ 0.5555555555555556 0.8888888888888889 0.5555555555555556 ]);
    elseif (order==4)
        param_coords = vec([ -0.86113631159405  -0.33998104358486   0.33998104358486   0.86113631159405]);
        weights  = vec([ 0.34785484513745   0.65214515486255   0.65214515486255   0.34785484513745]);
    elseif (order==5)
        param_coords = vec([ -0.906179845938664        -0.538469310105683         0.000000000000000         0.538469310105683         0.906179845938664])
        weights  = vec([0.236926885056189        0.478628670499367        0.568888888888889        0.478628670499367        0.236926885056189])
    else
        param_coords, weights  = opengauss(order)
    end
    return length(param_coords), param_coords, weights
end


function _triangle(npts=1)
    if npts == 1 # integrates exactly linear polynomials
        param_coords = [1.0/3. 1.0/3.];
        weights = reshape([1.0]/2.0,1,1);
    elseif npts == 3 # integrates exactly quadratic polynomials
        param_coords = [ 2.0/3 1.0/6; 1.0/6 2.0/3; 1.0/6 1.0/6 ];
        weights = [1.0/3 1.0/3 1.0/3]/2;
    elseif npts == 4 # integrates exactly quadratic polynomials
        param_coords = [0.333333333333333   0.333333333333333
                        0.200000000000000   0.200000000000000
                        0.600000000000000   0.200000000000000
                        0.200000000000000   0.600000000000000];
        weights = [-0.281250000000000
                    0.260416666666667
                    0.260416666666667
                    0.260416666666667];
    elseif npts == 6 # integrates exactly quartic polynomials
        param_coords = [ 0.816847572980459 0.091576213509771;
                        0.091576213509771 0.816847572980459;
                        0.091576213509771 0.091576213509771;
                        0.108103018168070 0.445948490915965;
                        0.445948490915965 0.108103018168070;
                        0.445948490915965 0.445948490915965];
        weights = [0.109951743655322*[1, 1, 1] 0.223381589678011*[1, 1, 1] ]/2;
    elseif npts == 7 # integrates exactly ? polynomials
        param_coords = [0.101286507323456   0.101286507323456
                        0.797426958353087   0.101286507323456
                        0.101286507323456   0.797426958353087
                        0.470142064105115   0.470142064105115
                        0.059715871789770   0.470142064105115
                        0.470142064105115   0.059715871789770
                        0.333333333333333   0.333333333333333];
        weights = [0.062969590272414
                    0.062969590272414
                    0.062969590272414
                    0.066197076394253
                    0.066197076394253
                    0.066197076394253
                    0.112500000000000];
    elseif npts == 9 # integrates exactly ? polynomials
        param_coords = [   0.437525248383384   0.437525248383384
            0.124949503233232   0.437525248383384
            0.437525248383384   0.124949503233232
            0.165409927389841   0.037477420750088
            0.037477420750088   0.165409927389841
            0.797112651860071   0.165409927389841
            0.165409927389841   0.797112651860071
            0.037477420750088   0.797112651860071
            0.797112651860071   0.037477420750088];
        weights = [   0.205950504760887
            0.205950504760887
            0.205950504760887
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223
            0.063691414286223] ./ 2;
    elseif npts == 12 # integrates exactly ? polynomials
        param_coords = [   0.063089014491502   0.063089014491502
            0.873821971016996   0.063089014491502
            0.063089014491502   0.873821971016996
            0.249286745170910   0.249286745170910
            0.501426509658179   0.249286745170910
            0.249286745170910   0.501426509658179
            0.310352451033785   0.053145049844816
            0.053145049844816   0.310352451033785
            0.636502499121399   0.310352451033785
            0.310352451033785   0.636502499121399
            0.053145049844816   0.636502499121399
            0.636502499121399   0.053145049844816];
        weights = [0.050844906370207
            0.050844906370207
            0.050844906370207
            0.116786275726379
            0.116786275726379
            0.116786275726379
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374
            0.082851075618374] ./ 2;
    elseif npts == 13 # integrates exactly ? polynomials
        param_coords = [0.333333333333333  0.333333333333333
                        0.479308067841923  0.260345966079038
                        0.260345966079038  0.479308067841923
                        0.260345966079038  0.260345966079038
                        0.869739794195568  0.065130102902216
                        0.065130102902216  0.869739794195568
                        0.065130102902216  0.065130102902216
                        0.638444188569809  0.312865496004875
                        0.638444188569809  0.048690315425316
                        0.312865496004875  0.638444188569809
                        0.312865496004875  0.048690315425316
                        0.048690315425316  0.638444188569809
                        0.048690315425316  0.312865496004875
        ];
        weights = [ -0.149570044467670
                    0.175615257433204
                    0.175615257433204
                    0.175615257433204
                    0.053347235608839
                    0.053347235608839
                    0.053347235608839
                    0.077113760890257
                    0.077113760890257
                    0.077113760890257
                    0.077113760890257
                    0.077113760890257
                    0.077113760890257
        ]'/2;
    else
        #nothing doing: this input is wrong
        error( "Unknown number of integration points $(npts)" )
    end
    return npts, param_coords, weights
end

function _tetrahedron(npts=1)
    if npts == 1 # integrates exactly linear polynomials
        param_coords = reshape([0.25,0.25,0.25],1,3);
        weights = reshape([1.0]/6.0,1,1);
    elseif npts == 4 # integrates exactly quadratic polynomials
        param_coords = [[0.13819660 0.13819660 0.13819660];
        [0.58541020 0.13819660 0.13819660];
        [0.13819660 0.58541020 0.13819660];
        [0.13819660 0.13819660 0.58541020]];;
        weights = [ 0.041666666666666666667   0.041666666666666666667   0.041666666666666666667   0.041666666666666666667];
    elseif npts == 5 #  Zienkiewicz #3.
        a =   1.0 / 6.0;
        b =   0.25;
        c =   0.5;
        d = - 0.8;
        e =   0.45;
        param_coords = [[b b b];
        [c a a];
        [a c a];
        [a a c];
        [a a a]];
        weights = [d  e  e  e  e]/6;
    else
        #nothing doing: this input is wrong
        error( "Unknown number of integration points $(npts)" )
    end
    return npts, param_coords, weights
end

"""
    quadrature(::Type{RefShapeInterval}, quadraturesettings = (kind = :default,))

Create a quadrature rule for the reference shape of an interval.

The default is Gauss integration rule, where the order is set with the keyword
`order`.
"""
function quadrature(::Type{RefShapeInterval}, quadraturesettings = (kind = :default,))
	kind = :default
	for apair in pairs(quadraturesettings)
	    sy, val = apair
	    if sy == :kind
	        kind = val
	    end
	end
    if (kind == :Gauss) || (kind == :default)
        # Extract arguments
        order = 1; 
        for apair in pairs(quadraturesettings)
            sy, val = apair
            if sy == :order
                order = val
            end
        end
        npts, param_coords, weights = _gauss1(order)
        return IntegRule(npts, reshape(param_coords, size(param_coords, 1), 1), vec(weights))
    else
        error("Integration rule $(kind) not available")
    end
end

"""
    quadrature(::Type{RefShapeTriangle}, quadraturesettings = (kind = :default,))

Create a quadrature rule for the reference shape of an triangle.

The default is a triangle rule, distinguished by the number of points set with
the keyword `npts`.
"""
function quadrature(::Type{RefShapeTriangle}, quadraturesettings = (kind = :default,))
	kind = :default
	for apair in pairs(quadraturesettings)
	    sy, val = apair
	    if sy == :kind
	        kind = val
	    end
	end
    if kind == :default
        # Extract arguments
        npts = 1; 
        for apair in pairs(quadraturesettings)
            sy, val = apair
            if sy == :npts
                npts = val
            end
        end
        npts, param_coords, weights = _triangle(npts)
        return IntegRule(npts, reshape(param_coords, size(param_coords, 1), 2), vec(weights))
    else
        error("Integration rule $(kind) not available")
    end
end

"""
    quadrature(::Type{RefShapeSquare}, quadraturesettings = (kind = :default,))

Create a quadrature rule for the reference shape of a square.

The default is Gauss integration rule, where the order is set with the keyword
`order`.
"""
function quadrature(::Type{RefShapeSquare}, quadraturesettings = (kind = :default,))
	kind = :default
	for apair in pairs(quadraturesettings)
	    sy, val = apair
	    if sy == :kind
	        kind = val
	    end
	end
    if (kind == :Gauss) || (kind == :default)
        # Extract arguments
        order = 1; 
        for apair in pairs(quadraturesettings)
            sy, val = apair
            if sy == :order
                order = val
            end
        end
        np, pc, w = _gauss1(order)
        param_coords = zeros(eltype(w),np^2,2);
        weights = zeros(eltype(w),np^2);
        r=1
        for i in 1:order
        	for j in 1:order
        		param_coords[r,:] = [pc[i] pc[j]];
        		weights[r] = w[i]*w[j];
        		r=r+1
        	end
        end
        npts = r - 1 # back off one to reach the total number of points
        return IntegRule(npts, param_coords, weights)
    else
        error("Integration rule $(kind) not available")
    end
end

"""
    quadrature(::Type{RefShapeTetrahedron}, quadraturesettings = (kind = :default,))

Create a quadrature rule for the reference shape of a tetrahedron.

The default is a one-point tetrahedron rule; other rules may be chosen based on
the number of points set with the keyword `npts`.
"""
function quadrature(::Type{RefShapeTetrahedron}, quadraturesettings = (kind = :default,))
    kind = :default
    for apair in pairs(quadraturesettings)
        sy, val = apair
        if sy == :kind
            kind = val
        end
    end
    if kind == :default
        # Extract arguments
        npts = 1; 
        for apair in pairs(quadraturesettings)
            sy, val = apair
            if sy == :npts
                npts = val
            end
        end
        npts, param_coords, weights = _tetrahedron(npts)
        return IntegRule(npts, reshape(param_coords, size(param_coords, 1), 3), vec(weights))
    else
        error("Integration rule $(kind) not available")
    end
end

end # module
