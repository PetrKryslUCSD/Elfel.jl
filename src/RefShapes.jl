module RefShapes


abstract type AbstractRefShape{MANIFDIM} end

struct RefShapePoint <: AbstractRefShape{0} end
struct RefShapeInterval <: AbstractRefShape{1} end
struct RefShapeSquare <: AbstractRefShape{2} end
struct RefShapeCube <: AbstractRefShape{3} end
struct RefShapeTriangle <: AbstractRefShape{2} end
struct RefShapeTetrahedron <: AbstractRefShape{3} end

"""
    manifdim(rs)

Get the manifold dimension of the reference shape.
"""
manifdim(::Type{T}) where {T<:AbstractRefShape{MANIFDIM}} where {MANIFDIM} = MANIFDIM
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
    elseif (order==6)
        param_coords = vec([-0.932469514203152
            -0.661209386466264
            -0.238619186083197
            0.238619186083197
            0.661209386466264
            0.932469514203152]);
        weights  = vec([0.171324492379171
            0.360761573048139
            0.467913934572691
            0.467913934572692
            0.360761573048139
            0.171324492379170]);
    elseif (order==7)
        param_coords = vec([-0.949107912342758
            -0.741531185599394
            -0.405845151377397
            -0.000000000000000
            0.405845151377397
            0.741531185599395
            0.949107912342758]);
        weights  = vec([0.129484966168870
            0.279705391489277
            0.381830050505119
            0.417959183673469
            0.381830050505119
            0.279705391489276
            0.129484966168870]);
    elseif (order==8)
        param_coords = vec([-0.960289856497536
            -0.796666477413627
            -0.525532409916329
            -0.183434642495650
            0.183434642495650
            0.525532409916329
            0.796666477413627
            0.960289856497536]);
        weights  = vec([0.101228536290376
            0.222381034453374
            0.313706645877887
            0.362683783378362
            0.362683783378362
            0.313706645877887
            0.222381034453374
            0.101228536290376]);
    elseif (order==9)
        param_coords = vec([-0.968160239507626
            -0.836031107326636
            -0.613371432700590
            -0.324253423403809
            0.000000000000000
            0.324253423403809
            0.613371432700591
            0.836031107326636
            0.968160239507626]);
        weights  = vec([0.081274388361574
            0.180648160694858
            0.260610696402935
            0.312347077040002
            0.330239355001259
            0.312347077040002
            0.260610696402936
            0.180648160694857
            0.081274388361574]);
    elseif (order==10)
        param_coords = vec([-0.973906528517171
            -0.865063366688985
            -0.679409568299025
            -0.433395394129247
            -0.148874338981631
            0.148874338981631
            0.433395394129247
            0.679409568299024
            0.865063366688984
            0.973906528517172]);
        weights  = vec([0.066671344308688
            0.149451349150581
            0.219086362515981
            0.269266719309996
            0.295524224714752
            0.295524224714753
            0.269266719309996
            0.219086362515982
            0.149451349150581
            0.066671344308688]);
    else
        error("Gauss rule of order $(order) not available")
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

end # module
