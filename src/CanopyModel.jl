

###############################################################################
############################## Canopy structure ###############################
###############################################################################

# Store all the vectors and matrices associated to the canopy
struct Canopy{N, N1, N12}
    ΔL::Float64
    Lₘ::SVector{N,Float64}
    Lᵤ::SVector{N,Float64}
    Fu::SMatrix{N1, N1, Float64, N12}
    Fd::SMatrix{N1, N1, Float64, N12}
    fdif::SVector{N, Float64}
    Vcmax25::SVector{N, Float64}
    Rd25::SVector{N, Float64}
    Jmax25::SVector{N, Float64}
    α::SVector{N, Float64}
    sigma::SVector{N1, Float64}
end

# Create the canopy structure for each day (pre-allocations and pre-calculations)
function Canopy(pars, vars, ::Val{N} = Val(15)) where N
    # Leaf area index per layer
    ΔL = vars.GLAI/N
    # Intermediate variables to be reused
    Lₘ   = SVector{N, Float64}((i - 1)*ΔL + ΔL/2 for i in 1:N)
    Lᵤ   = SVector{N, Float64}(i*ΔL for i in 1:N)
    # Update sky model with the extinction coefficients for diffuse light
    updated_sky = ksky(pars.angles, pars.sky) 
    fdif = calc_fdiff(Lᵤ, updated_sky)
    # Matrices for redistribution of scattered radiation
    nθ = length(updated_sky.β)
    Fu, Fd = calc_fscat(ΔL,nθ, Lᵤ, pars.angles)
    # Distribution of total nitrogen and N partitioning coefficients within canopy
    Np = @. 0.65*(vars.Nlmax*exp(-vars.kN*Lₘ) - vars.Nmin)
    f_Nc = @. vars.f_Ncmn + (vars.f_Ncmx - vars.f_Ncmn)*(Lₘ/vars.GLAI)^vars.pf_Nc
    f_Nr = @. vars.f_Nrmx - (vars.f_Nrmx - vars.f_Nrmn)*(Lₘ/vars.GLAI)^vars.pf_Nr
    f_No = 1.0 .- f_Nc .- f_Nr
    f_Nt = f_No./(1.0 .+ pars.aj.*pars.ks)
    # Compute photosynthetic traits
    Nc = f_Nc.*Np
    α = Nc./(Nc .+ pars.Knc)
    Vcmax25 = pars.ar.*f_Nr.*Np
    Jmax25 = pars.aj.*f_Nt.*Np
    Rd25 = pars.f_Rd.*Vcmax25
    # Scattering coefficient of each layer from chlorophyll content (+ soil)
    sigma = vcat(1.0 .- α, pars.sigma_soil)
    # Construct the object
    Canopy{N, N + 1, (N + 1)*(N + 1)}(ΔL, Lₘ, Lᵤ, Fu, Fd, fdif, Vcmax25,  Rd25, Jmax25, α, sigma)
end

# Extinction coefficients to different sectors of the sky given a leaf angle distribution
# and Gaussian-Legendre integrator (not needed for some distributions). This function also
# returns the values contained inside sky to be reused later on. 
function ksky(angles, sky) 
    ks = SVector{length(sky.β), Float64}(k(β, angles) for β in sky.β)
    (k = ks, β = sky.β, f = sky.f)    
end


#### Precalculations for scattered light ####

"""
calc_fscat(ΔL, nL, sky)
Matrix with fraction of scattered irradiance emmitted by a layer that reaches another layer
calculated using the modified angular distribution of scattered radiation (Goudriaan, 1977, eq. 2.35)

# Examples
```julia
```
"""
function calc_fscat(ΔL, nθ, Lᵤ::SVector{N, Float64}, angles) where N
    @inbounds begin
        # Angular distribution for scattered light (modified by leaf angle)
        scat = scattered_model(angles, nθ, ΔL[1])
        scat_k = ksky(angles, scat) 
        # Fraction of scattered radiation that reaches a particular distance in leaf area index
        fI = SVector{N + 1, Float64}(1.0, (sum(@. exp(-scat_k.k*L)*scat_k.f) for L in Lᵤ)...)
        # Use fI to compute the fraction of scattered radiation intercepted by each layer
        fdiff = fI[1:end-1] .- fI[2:end]
        # Arrange the fractions into a matrix for easier computation
        F = @MMatrix zeros(N + 1, N + 1)
        for i in 2:N
            F[:,i] = vcat(fdiff[i-1:-1:1], 0.0, fdiff[1:end-i+1])
        end
        F[:,1] = vcat(0.0, fdiff)
        F[:,end] = vcat(fI[end-1:-1:1], 0.0)
        F = transpose(F)
        Fu = SMatrix{N + 1, N + 1, Float64}(UpperTriangular(F))
        Fd = SMatrix{N + 1, N + 1, Float64}(LowerTriangular(F))
        return (Fu = Fu, Fd = Fd)
    end
end

#### Precalculations for diffuse light ####

# Diffuse irradiance (relative to value on top of canopy) that is intercepted by each
# layer of the canopy
function calc_fdiff(Lᵤ::SVector{N, Float64}, sky) where N
    # Fraction of diffuse radiation that reaches a particular distance in leaf area index
    fI = SVector{N + 1, Float64}(1.0, (sum(@. exp(-sky.k*L)*sky.f) for L in Lᵤ)...)
    # Use fI to compute the fraction of diffuse intercepted by each layer
    @inbounds fI[1:end-1] .- fI[2:end]
end


###############################################################################
################################ Black leaves #################################
###############################################################################

# Compute average irradiance for each layer + irradiance reaching the soil surface
function black_canopy(can::Canopy{N, N1, N12}, env, pars, vars) where {N, N1, N12}
    # Calculate the diffuse light intercepted by each layer
    Id = env.Id0.*can.fdif
    # Calculate the average irradiance increase due to direct solar radiation
    kb    = calc_kdir(env.β, pars.angles)
    fsun  = exp.(.-kb.*can.Lₘ)
    f_dir = calc_fdir(can, kb)
    Ib    = env.Ib0.*f_dir
    # Calculate irradiance reaching the soil surface
    @inbounds Idsoil  = env.Id0*sum(exp(-pars.sky.k[i]*vars.GLAI)*pars.sky.f[i] 
                                    for i in 1:length(pars.sky.f))
    Ibsoil  = env.Ib0*exp(-kb*vars.GLAI)
    Isoil   = Idsoil + Ibsoil
    # Average irradiance intercepted by leaves in each layer
    @inbounds Ip = SVector{N1, Float64}((Id[i] + Ib[i] for i in 1:N)..., Isoil)

    return Ip, Id, Ib, fsun
end

# Compute the extinction coefficient to direct solar radiation with angle β
calc_kdir(β, angles) = k(β, angles)

# Compute the fraction of direct solar radiation on top of canopy that is 
# absorbed by a particular canopy layer
function calc_fdir(can::Canopy{N, N1, N12}, kb) where {N, N1, N12}
    fI = SVector{N + 1, Float64}(1.0, (exp(-kb*L) for L in can.Lᵤ)...)
    @inbounds SVector{N, Float64}(fI[1:N] .- fI[2:N1])
end


###############################################################################
################################# Grey leaves #################################
###############################################################################

# Update incoming scattered irradiance for each layer
function update_scatter(can::Canopy{N,N1, N12}, Sio, SIp)::Tuple{SVector{N1, Float64}, Float64} where {N, N1, N12}
    # Compute scattered radiation emitted by each layer (from previous iteration)
    So  = SIp .+ can.sigma.*Sio./2.0
    # New intercepted scattered radiation by each layer
    Siu = can.Fu*So
    Sid = can.Fd*So
    Si  = Siu .+ Sid
    # Error norm for following convergence
    n = sum(abs.((Si .- Sio)./Si))
    return Si, n
end


# Iteratively compute scattered irradiance
function calc_scatter(can::Canopy{N,N1, N12}, Ip)::Tuple{SVector{N, Float64}, Float64} where {N,N1, N12}
    # Scattering from incident primary radiation
    SIp = Ip.*can.sigma./2.0
    # Resolve secondary scattering iteratively
    Si =  @SVector zeros(N + 1)
    Si, n = update_scatter(can, Si, SIp)
    for _ in 1:50
        Si, n = update_scatter(can, Si, SIp)
        n < 1e-4 && break
    end
    return @inbounds Si[1:N], Si[end]
end


# Update incident irradiance on the leaves by scattering
function grey_canopy(can::Canopy{N,N1, N12}, env, pars, vars) where {N,N1, N12}
    # Compute intercepted radiation assuming black leaves
    Ip, Id, Ib, fsun = black_canopy(can, env, pars, vars)
    # Calculate additional intercepted radiation due to scattering (includes soil reflection)
    Si, Sisoil = calc_scatter(can, Ip)
    Ishade = Id .+ Si
    Ilayer = Ishade .+ Ib
    return Ilayer, Ishade, Ib, fsun, Ip[end] + Sisoil
end

###############################################################################
############################### Photosynthesis ################################
###############################################################################

# Compute the projection of beam solar irradiance onto a leaf angle
function kbeam(λ, β)
    if β >= λ  
        cos(λ) 
    else     
        2/π*(cos(λ)*asin(tan(β)/tan(λ)) + sqrt(sin(λ)^2 - sin(β)^2)/sin(β))
    end
end

# Compute CO2 assimilation of a sunlit leaf based on the relevant solar angles
function A_sun_ang(ang, β, Omega, Ib, Ishade, pars, Rd, k2ll, Jmax, gamma_star, gm, fvpd, Ca)::Float64
    # The angles of the leaf
    @inbounds λ = ang[1]
    @inbounds α = ang[2]
    # The irradiance incident on the leaf
    Iinc = Ishade + Ib*t(β, Omega, λ, α)
    # Compute gross CO2 assimilation with the model for the incident irradiance
    x1_j = J(k2ll, pars.theta, Jmax, Iinc)
    Aj = CalcAnC3(gm, pars.gs0, fvpd, pars.gb, 2.0*gamma_star, x1_j, gamma_star, Rd, Ca)
    # Multiply by the probability density of (λ, α) under the leaf angle distribution
    Aj*pars.angles.fλ(λ)/2/π
end

# Integrate photosynthesis in the sunlit fraction over the leaf angle distribution
function calc_Asun(env, pars, Ib, Ishade, Ac, Rd, k2ll, Jmax, gamma_star, gm, fvpd)
    f = ang -> A_sun_ang(ang, env.β, env.Omega, Ib, Ishade, pars, Rd, k2ll, Jmax, gamma_star, gm, fvpd, env.Ca)
    Aj = hcubature(f, SVector(0.0, 0.0), SVector(π/2, 2π), rtol = pars.angles.rtol, atol = 0.0, maxevals = 5_000)[1]
    min(Ac, Aj) + Rd
end

# Given irradiances and photosynthetic parametes per layer, compute CO2 assimilation for shaded and sunlit fractions
function layer_assimilation(can::Canopy{N,N1, N12}, env, pars, Ishade, Ib) where {N,N1, N12}

    k2ll, Jmax, Vcmax, Kmapp, Rd, gamma_star, gm, fvpd = 
            temperature_correction(can.α, can.Jmax25, can.Rd25, can.Vcmax25, env.Tleaf, env.VPD, pars)

    # Photosynthesis limited by Rubisco per layer
    Ac = CalcAnC3.(gm, pars.gs0, fvpd, pars.gb, Kmapp, Vcmax, gamma_star, Rd, env.Ca)

    # Calculate photosynthesis of shaded leaves
    x1_j = J.(k2ll, pars.theta, Jmax, Ishade)
    Aj = CalcAnC3.(gm, pars.gs0, fvpd, pars.gb, 2.0*gamma_star, x1_j, gamma_star, Rd, env.Ca)
    Ashade  = min.(Ac, Aj) .+ Rd

    # Calculate photosynthesis of sunlit leaves
    @inbounds Asun = SVector{N, Float64}(calc_Asun(env, pars, Ib, Ishade[i], Ac[i], Rd[i], k2ll[i], Jmax[i],
                                                   gamma_star, gm, fvpd) for i in 1:N)
    return Ashade, Asun
end

# Compute canopy photosynthesis
function Acan(can::Canopy{N,N1, N12}, env, pars, vars) where {N,N1, N12}
    # Compute irradiance intercepted by each layer
    Ilayer, Ishade, _, fsun, _ = grey_canopy(can, env, pars, vars)
    # Compute CO2 assimilation/leaf area in each layer by separating shaded and sunlit fractions
    Ashade, Asun = layer_assimilation(can, env, pars, Ishade, env.Ib0/sin(env.β))
    # Scale up assimilation by leaf area in each layer and fraction
    Lsun   = can.ΔL.*fsun
    Lshade = can.ΔL.*(1.0 .- fsun)
    Alayer = Ashade.*Lshade .+ Asun.*Lsun
    # Canopy CO2 assimilation
    return sum(Alayer), Alayer, Ilayer
end

# Compute total daily assimilation in mol/m2
function calc_Adiurnal(weather, can, pars,vars, method)
    f = t -> Acan(can, interpolate_meteo(weather, t), pars, vars)[1]
    Aday = integrate(method, f, 1e-6, weather.DL - 1e-6)
    # Convert from mol/m2/s to mol/m2
    return Aday*3600.0
end
