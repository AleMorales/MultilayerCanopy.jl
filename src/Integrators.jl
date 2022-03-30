# Objects to store nodes from numerical quadrature rules
abstract type Integrator end


#############################################################
############ Fixed Gaussian-Legendre integrator #############
#############################################################

# Fixed-grid gaussian-Legendre integrator
struct GLIntegrator{N} <: Integrator
  nodes::SVector{N,Float64}
  weights::SVector{N,Float64}
end

"""
  GLIntegrator(n)
Create object to store nodes and weights (w) of a Gaussian-Legendre integrator of
order n, normalize to the interval [0,1] and with ∑w = 1

# Examples
```julia
GL = GLIntegrator(5)
GL.nodes
GL.weights
```
"""
function GLIntegrator(::Val{N}) where N
  nodes, weights = gausslegendre(N)
  nodes .= (nodes .+ 1.0)./2.0
  weights ./= 2.0
  GLIntegrator(SVector{N,Float64}(nodes), SVector{N,Float64}(weights))
end
Base.size(g::GLIntegrator{N}) where N = N

function integrate(method::GLIntegrator{N}, f, a, b) where {N}
    out = 0.0
    @inbounds for i in 1:N
        out += f(a + method.nodes[i]*(b - a))*method.weights[i]
    end
    out*(b - a)
end


#########################################################################
############ Fixed Partition + Gaussian-Legendre integrator #############
#########################################################################

"""
  PGLIntegrator(Val(n), Val(p))
Create object to store nodes and weights (w) of a Gaussian-Legendre integrator of
order n to be evaluated on p sub-divisions of the original range.

# Examples
```julia
GL = PGLIntegrator(Val(5), Val(3))
GL.nodes
GL.weights
```
"""
function PGLIntegrator(::Val{N}, ::Val{P}) where {N, P}
  # Original N-order grid
  nodes_sub, weights_sub = gausslegendre(N)
  # Scale down the weights and repeat as needed
  weights_sub ./= (2P)
  weights = SVector{N*P, Float64}(repeat(weights_sub, P)...)
  # Normalize the nodes
  nodes_sub .= (nodes_sub .+ 1.0)./2.0
  lower  = repeat([(i - 1)/P for i in 1:P], inner = N) 
  nodes  = repeat(nodes_sub, P)./P .+ lower
  # Create the integrator object with the right nodes and weights
  GLIntegrator(SVector{N*P,Float64}(nodes),weights)
end


#############################################################
########### Adaptive Gaussian-Legendre integrator ###########
#############################################################

struct AdaptGL <: Integrator 
    rtol::Float64
end

function integrate(method::AdaptGL, f, a, b)
    quadgk(f, a, b, rtol = method.rtol)[1]
end

#############################################################
################ Fixed midpoint integrator ##################
#############################################################

# Fixed-grid gaussian-Legendre integrator
struct MidPoint{N} <: Integrator end
MidPoint(n) = MidPoint{n}()

function integrate(method::MidPoint{N}, f, a, b) where {N}
    sum(f((b - a)*(2i - 1)/(2N) + a) for i in 1:N)*(b - a)/N
end

#############################################################
################ Fixed trapezoid integrator #################
#############################################################

struct Trapezoid{N} <: Integrator end
Trapezoid(n) = Trapezoid{n}()

function integrate(method::Trapezoid{N}, f, a, b) where {N}
    (sum(f((b-a)*i/N + a) for i in 1:N-1) +
     (f(a) + f(b))/2)*(b - a)/N
end
