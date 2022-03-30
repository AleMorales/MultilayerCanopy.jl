
"""
    LeafAngleModel(fλ, rtol)

Combination of a leaf angle distribution object (fλ) and the relative tolerance
to be used by the different numerical algorithms that may be used (e.g. to compute
extinction coefficients or average sunlit photosynthesis.)

# Examples
```julia
lm = LeafAngleModel(Ellipsoidal(1), 1e-4)
```
"""
struct LeafAngleModel{L}
    fλ::L
    rtol::Float64
    LeafAngleModel(fλ::L, rtol = 1e-4) where L = new{L}(fλ, rtol)
end

####################################
##### Leaf angle distributions #####
####################################

abstract type LeafAngle end

"""
    Planophile()
Functor that implements the planophile leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
p = Planophile()
p(pi/2)
```
"""
struct Planophile <: LeafAngle end
(p::Planophile)(λ) = 2/π*(1 + cos(2*λ))

"""
    Erectophile()
Functor that implements the erectophile leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
e = Erectophile()
e(pi/2)
```
"""
struct Erectophile <: LeafAngle end
(p::Erectophile)(λ) = 2/π*(1 - cos(2*λ))

"""
    Plagiophile()
Functor that implements the plagiophile leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
p = Plagiophile()
p(pi/2)
```
"""
struct Plagiophile <: LeafAngle end
(p::Plagiophile)(λ) = 2/π*(1 - cos(4*λ))

"""
    Extremophile()
Functor that implements the extremophile leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
e = Extremophile()
e(pi/2)
```
"""
struct Extremophile <: LeafAngle end
(p::Extremophile)(λ) = 2/π*(1 + cos(4*λ))

"""
    Uniform()
Functor that implements the uniform leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
u = Uniform()
u(pi/2)
```
"""
struct Uniform <: LeafAngle end
(p::Uniform)(λ) = 2/π

"""
    Spherical()
Functor that implements the spherical leaf angle distribution from Table 2
in Weiss et al. (2004)

# Examples
```julia
s = Spherical()
s(pi/2)
```
"""
struct Spherical <: LeafAngle end
(p::Spherical)(λ) = sin(λ)

"""
    Beta(μ, ν)
Functor that implements the beta leaf angle distribution from Wang et al. (2007)
where μ and ν are two empirical parameters

# Examples
```julia
b = Beta(2.770, 1.172)
b(pi/2)
```
"""
struct Beta <: LeafAngle
  μ::Float64
  ν::Float64
end
(b::Beta)(λ) = (t = 2/pi*λ; ((1 - t)^(b.μ - 1)*t^(b.ν - 1))/beta(b.μ, b.ν))

"""
    Ellipsoidal(χ)
Functor that implements the ellipsoidal leaf angle distribution from Campbell (1986)
where χ is an empirical parameter

# Examples
```julia
e = Ellipsoidal(1)
e(pi/2)
```
"""
struct Ellipsoidal <: LeafAngle
  χ::Float64
  Λ::Float64
end
function Ellipsoidal(χ)
  if χ < 1
    ϵ = sqrt(1 - χ*χ)
    Λ = χ + asin(ϵ)/ϵ
  elseif χ == 1
    Λ = 2
  else
    ϵ = sqrt(1 - 1/(χ*χ))
    Λ = χ + log((1 + ϵ)/(1 - ϵ))/(2*ϵ*χ)
  end
  Ellipsoidal(χ, Λ)
end
function (e::Ellipsoidal)(λ)
  sinλ = sin(λ)
  cosλ = cos(λ)
  den = (cosλ*cosλ + e.χ*e.χ*sinλ*sinλ)
  2*e.χ*e.χ*e.χ*sinλ/(e.Λ*den*den)
end

####################################
###### Extinction coefficients #####
####################################

"""
    G(β, λ)
Projection of leaves with inclination angle λ on direction with inclination angle β
averaged over all leaf azimuth angles assuming uniform leaf azimuth distribution.
Eqn 2.3 from Goudriaan (1977) or Eqn 8 from Wang et al. (2007).

# Examples
```julia
G(pi/2, pi/2)
G(pi/2, 0)
```
"""
function G(β, λ)
  if λ ≤ β
    sin(β)cos(λ)
  else
    sinβ = sin(β)
    cosβ = cos(β)
    tanβ = sinβ/cosβ
    sinλ = sin(λ)
    cosλ = cos(λ)
    tanλ = sinλ/cosλ
    2/π*(sinβ*cosλ*asin(tanβ/tanλ) + sqrt(sinλ*sinλ - sinβ*sinβ))
  end
end


"""
    k(β, L)
Extinction coefficient for inclination angle β by averaging over all leaf
azimuth and inclination angles assuming uniform leaf azimuth distribution and fλ
inclination angle distribution.

Implements Eqn 4 & 7 from Wang et al (2007) or Eqn 2.4 & 2.28 from Goudriaan (1977)
but using adaptive Gaussian-Kronrod integration (the relative tolerance is set by
the corresponding argument in the `LeafAngleModel` object).

For Spherical and Ellipsoidal distributions, analytical solutions from Goudriaan
(1977) and Campbell(1986) are used, respectively.

# Examples
```julia
# Spherical distribution
L = LeafAngleModel(Spherical())
β = pi/4
ks = k(β, L)
# Ellipsoidal distribution with χ = 1
L = LeafAngleModel(Ellipsoidal(1.0))
ke = k(β, L)
# Beta distribution approximating the Spherical distribution
t̄ = 0.6371907
vt = 0.1600465
μ = (1 - t̄)*(t̄*(1 - t̄)/vt - 1)
ν = t̄/(1 - t̄)*μ
L = LeafAngleModel(Beta(μ, ν))
kb = k(β, L)
```
"""
function k(β, L::LeafAngleModel)
  out = quadgk(λ -> G(β, λ)*L.fλ(λ), 1e-6, π/2 - 1e-6, rtol = L.rtol)[1]
  out*π/2.0/sin(β)
end

# Analytical solutions for spherical and Ellipsoidal distributions
k(β, L::LeafAngleModel{Spherical}) = 0.5/sin(β)
function k(β, L::LeafAngleModel{Ellipsoidal})
  tanβ = tan(β)
  sqrt(L.fλ.χ*L.fλ.χ + 1/(tanβ*tanβ))/L.fλ.Λ
end


####################################
######### Other expressions ########
####################################

"""
    scattered_model(L, nθ, ΔL)
Calculate angular distribution of scattered radiation as a modification of the 
uniform distribution, where L is the `LeafAngleModel` object, nθ is the number of
angles in which we will subdivide the hemisphere and ΔL is the leaf area of a layer
in the canopy. This is an implementation of Eqn. 2.35 by Goudriaan (1977).
 
# Examples
```julia
L = LeafAngleModel(Spherical())
scattered_model(L, 10, 0.1)
```
"""
function scattered_model(L::LeafAngleModel, nθ, ΔL)
    Δθ = π/2/nθ
    θₗ = 0.0:Δθ:π/2 - Δθ
    θᵤ = Δθ:Δθ:π/2
    βs = @. π/2 - (θₗ + θᵤ)/2
    ks = [k(β, L) for β in βs]
    f = @. (cos(θₗ)^2 - cos(θᵤ)^2)*exp(-ks*ΔL)
    f /= sum(f)
    return (β = βs, f = f)
end


"""
    t(β, ϕ, λ, α)
Cosine of the angle of incidence between a solar ray with elevation β and 
azimuth ϕ and a leaf with elevation angle λ and azimuth orientation α. Notice
that we take the absolute value as the sign of the expression depends on the
side of the leaf (we assume same rate of photosynthesis on both sides). This
implements Eqn. 2.1 by Goudriaan (1977).

# Examples
```julia
t(π/4, π/3, π/5, 0.0)
```
"""
t(β, ϕ, λ, α) = abs(sin(β)*cos(λ) + cos(β)*sin(λ)*cos(α - ϕ))

