#### Misc Data Types ####
# const ThreeVec{T} = SVector{3,T}

#### Misc Aux Functions ####

## Range equally spaced in log10 scale ##
function GeomRange(start, stop, len::Integer)
  return 10 .^ (LinRange(log10(start), log10(stop), len))
end

## x*coth(x) function with special treatment of x=0 ##
function xcoth(x)
  if x == 0.0
    return 1.0
  else
    return x*coth(x)
  end
end

## Sinh(x)/x ##
sinhc(x::Real) = real(sinc(1im*x/π))

## Coth(x) - 1/x ##
cothminv(x::Number) = iszero(x) ? zero(x) : coth(x) - one(x)/x

## Coth(x)/x - 1/x^2 ##
cothinvminv2(x::Number) = iszero(x) ? one(x)/3 : (coth(x) - one(x)/x)/x

## Double integral with QuadGK ##
function dblquadgk(f, a::AbstractArray{T}, b::AbstractArray{T};
                   rtol=sqrt(eps(T)), atol=zero(T), maxevals=10^7, order=7) where T<:AbstractFloat
  J(x) = quadgk(y -> f(x,y), a[2], b[2], atol=atol, maxevals=maxevals, order=order)[1]
  K = quadgk(x -> J(x), a[1], b[1], atol=atol, maxevals=maxevals, order=order)[1]
  return K
end

function dblquadgk(f, a, b; kws...)
  T = promote_type(eltype(a), eltype(b))
  aT = [convert(T, a[1]), convert(T, a[2])]
  bT = [convert(T, b[1]), convert(T, b[2])]
  return dblquadgk(f, aT, bT; kws...)
end

## Integrate function over semi-real axis avoiding singularity ##
function integrate_semi_real(f, singularity = nothing)
  integrand_chvar = u -> f(u[1]/(1 - u[1]))/(1 - u[1])^2

  if isnothing(singularity)
    integrand_eval = quadgk(integrand_chvar, 0.0, 1.0)[1]
  else
    x0 = singularity[1]
    ε  = singularity[2]
    u0 = x0/(1 + x0)
    if isnan(u0)
      println(singularity)
    end
    integrand_eval_l = quadgk(integrand_chvar, 0.0, u0 - ε)[1]
    integrand_eval_r = quadgk(integrand_chvar, u0 + ε, 1.0)[1]
    integrand_eval   = integrand_eval_l + integrand_eval_r
  end

  return integrand_eval
end

## Spherical coordinates ##
@inline xsph(θ, φ) = sin(θ)*cos(φ)
@inline ysph(θ, φ) = sin(θ)*sin(φ)
@inline zsph(θ, φ) = cos(θ)
@inline dΩ(θ, φ) = sin(θ)
@inline Xsph(Ω) = [xsph(Ω...) ; ysph(Ω...) ; zsph(Ω...)]

@inline θsph(x, y, z) = acos(z)
@inline φsph(x, y, z) = iszero(x) && iszero(y) ? zero(x) : atan(y/x)
@inline θφsph(X) = [θsph(X...) ; φsph(X...)]
