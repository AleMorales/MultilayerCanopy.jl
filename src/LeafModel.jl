# Solution to system of equations involving CO2 sinks and diffusion from air to carboxylation site
function CalcAnC3(gm, gs0, fvpd, gb, x2, x1, gamma_star, Rd, Ca)::Float64
    a = gs0*(x2 + gamma_star) + (gs0/gm + fvpd)*(x1 - Rd)
    b = Ca*(x1 - Rd) - gamma_star*x1 - Rd*x2
    c = Ca + x2 + (1.0/gm + 1.0/gb)*(x1 - Rd)
    d = x2 + gamma_star + (x1 - Rd)/gm
    m = 1.0/gm + (gs0/gm + fvpd)*(1.0/gm + 1.0/gb)
    p = -(d + (x1 - Rd)/gm + a*(1.0/gm + 1.0/gb) + (gs0/gm + fvpd)*c)/m
    q = (d*(x1 - Rd) + a*c + (gs0/gm + fvpd)*b)/m
    r = (-a*b)/m
    U = (2.0*p^3.0 - 9.0*p*q + 27.0*r)/54.0
    Q = (p^2.0 - 3.0*q)/9.0
    psi = acos(U/sqrt(Q^3.0))
    A = -2.0*sqrt(Q)*cos(psi/3.0) - p/3.0
end

const R =  8310.0
const T0 =  298.15
const O2 = 210e-3
peaked(T, p25, E, S, D) = (p25*exp(((T - T0)*E)/(T0*R*T))*(1.0 + exp((T0*S - D)/(R*T0))))/(1.0 + exp((T*S - D)/(R*T)))
arrhenius(T, p25, E) = p25*exp(((T - T0)*E)/(T0*R*T))
parabolic(T, pOpt, Topt, Ω) = pOpt*exp(-(T - Topt)*(T - Topt)/(Ω*Ω))

# Apply tempeature correction
function temperature_correction(α, Jmax25, Rd25, Vcmax25, Tleaf, VPD, pars)
    # Effect of temperature of photosynthetic traits and other calculations
    Φ2 = parabolic(Tleaf, pars.Φ2LL, pars.Topt_Φ2, pars.Ω)
    s = Φ2*(1.0 - pars.fcyc)/(1.0 - pars.fcyc + Φ2/pars.Φ1LL)
    k2ll = α.*s
    Jmax = peaked.(Tleaf, Jmax25, pars.E_Jmax, pars.S_Jmax, pars.D_Jmax)
    Rd = arrhenius.(Tleaf, Rd25, pars.E_Rd)
    Vcmax = arrhenius.(Tleaf, Vcmax25, pars.E_Vcmax)
    Kmc   = arrhenius(Tleaf, pars.Kmc25, pars.E_Kmc)
    Kmo   = arrhenius(Tleaf, pars.Kmo25, pars.E_Kmo)
    Kmapp = Kmc*(1.0 + O2/Kmo)
    Sco   = arrhenius(Tleaf, pars.Sco25, pars.E_Sco)
    gamma_star = (0.5*O2)/Sco
    gm = peaked(Tleaf, pars.gm25, pars.E_gm, pars.S_gm, pars.D_gm)
    fvpd = 1.0/(1.0 + VPD/pars.D0)

    return k2ll, Jmax, Vcmax, Kmapp, Rd, gamma_star, gm, fvpd
end

# The hyperbol;ic light response curve
J(k, theta, Jmax, I) = (k*I + Jmax - sqrt((k*I + Jmax)^2.0 - 4.0*theta*k*Jmax*I))/(2.0*theta)/4

# INTERFACE: Compute net CO2 assimilation given parameter values and environmental conditions
# The canopy model does not use this function to avoid redundant computations
function Ag(pars, Np, f_Nc, f_Nr, VPD, PAR, Tleaf, Ca)

    # Compute the photosynthetic traits from nitrogen content and fractions
    f_No = 1.0 .- f_Nc .- f_Nr
    f_Nt = f_No./(1.0 .+ pars.aj.*pars.ks)
    Nc = f_Nc.*Np
    α = Nc./(Nc .+ pars.Knc)
    Vcmax25 = pars.ar.*f_Nr.*Np
    Jmax25 = pars.aj.*f_Nt.*Np
    Rd25 = pars.f_Rd.*Vcmax25

    # Apply temperature corrections and compute derived traits
    k2ll, Jmax, Vcmax, Kmapp, Rd, gamma_star, gm, fvpd = 
                        temperature_correction(α, Jmax25, Rd25, Vcmax25, Tleaf, VPD, pars)

    # Photosynthesis limited by Rubisco per layer
    Ac = CalcAnC3.(gm, pars.gs0, fvpd, pars.gb, Kmapp, Vcmax, gamma_star, Rd, Ca)

    # Photosynthesis limited by electron transport
    x1_j = J.(k2ll, pars.theta, Jmax, PAR)
    Aj = CalcAnC3.(gm, pars.gs0, fvpd, pars.gb, 2.0*gamma_star, x1_j, gamma_star, Rd, Ca)

    # Take the minimum of the two
    A = min.(Ac, Aj) .+ Rd

    return A, min.(Ac, Aj), Ac, Aj
end
