
# Parameters are constants determined by "crop genetics", general properties or settings
function parameters(;
    # Rubisco CO2/O2 specificity
    Sco25 = 2800.0, # Sc/o parameter at 25 C (bar/bar)
    E_Sco = -24.46, # Apparent activation energy of Sc/o (kJ/mol)
    # Michaelis-Menten constants Rubisco
    Kmc25 = 270.0, # Km for CO2 at 25 C (μmol/mol)
    E_Kmc = 80.99, # Activation energy of Kmc (kJ/mol)
    Kmo25 = 165.0, # Km for O2 at 25 C (mmol/mol)
    E_Kmo = 23.72, # Activation energy of Kmo (kJ/mol)
    # Rubisco activity
    ar =3.5/786, # Ratio between Vcmax25 and Nr (1/s)
    E_Vcmax = 65.33, # Activation energy of Vcmax (kJ/mol)
    # Electron transport
    theta = 0.7, # Curvature parameter
    Knc = 0.076/25.0, # Nc at which the leaf absorbs 50% of incoming PAR (mol/m2)
    Phi2LL = 0.85, # PSII quantum yield at low light
    Phi1LL = 1.0, # PSII quantum yield at low light
    fcyc = 0.1, # Fraction of cyclic electron transport
    Topt_Phi2 = 22.5 + 273.15, # Optimal temperature of PSII quantum yield
    Omega = 36.5, # sigma/√2 in the Gaussian temperature function for PSII quantum yield
    aj = 15.9e-3, # Ratio between Jmax25 and Nt (1/s)
    ks = 125.0, # Ratio between Ns and Jmax25 (s)
    E_Jmax = 30.0, # Activation energy Jmax (kJ/mol)
    D_Jmax = 200.0, # Deactivation energy of Jmax (kJ/mol)
    S_Jmax = 0.65, # Entropy coefficient of Jmax (kJ/mol/K)
    # Mitochondrial respiration
    E_Rd = 46.39, # Activation energy of Rd (kJ/mol)
    f_Rd  = 0.01, # Ratio between Rd25 and Vcmax25
    # Mesophyll conductance
    gm25 = 0.4, # Mesophyll conductance (mol/m2/s)
    E_gm = 70.2, # Activation energy gm (kJ/mol)
    S_gm = 0.32, # Entropy coefficient of gm (kJ/mol/K)
    D_gm = 94.0, # Deactivation energy of gm (kJ/mol)
    # Stomatal conductance
    gs0 = 0.05/1.56, # Minimum stomatal conductance (mol/m2/s)
    a1 = 0.9, # Empirical parameter capturing the effect of VPD on stomatal conductance
    b1 = 0.15, # Empirical parameter capturing the effect of VPD on stomatal conductance (1/kPa)
    # Boundary layer conductance
    gb = 0.5, # Boundary layer conductance (mol/m2/s)
    # Parameters related to canopy structure
    angles = LeafAngleModel(Spherical(), 1e-3),
    Ncmin   = 8.4e-3, # Minimum leaf nitrogen content (g N/g DW)
    # Other parameters related to light interception and climate
    sigma_soil = 0.21, # Albedo of the soil surface
    sky = ksky(angles, standard_sky(Val(10))),
    lat = 52.0,
    Ca = 4e2
)

(
    Sco25 = Sco25,
    E_Sco = E_Sco*1e3,
    Kmc25 = Kmc25*1e-6,
    E_Kmc = E_Kmc*1e3,
    Kmo25 = Kmo25*1e-3,
    E_Kmo = E_Kmo*1e3,
    ar = ar,
    E_Vcmax = E_Vcmax*1e3,
    theta = theta,
    Knc = Knc,
    Phi2LL = Phi2LL,
    Phi1LL = Phi1LL,
    fcyc = fcyc,
    Topt_Phi2 = Topt_Phi2,
    Omega = Omega,
    aj = aj,
    ks = ks,
    E_Jmax = E_Jmax*1e3,
    D_Jmax = D_Jmax*1e3,
    S_Jmax = S_Jmax*1e3,
    E_Rd = E_Rd*1e3,
    f_Rd  = f_Rd,
    gm25 = gm25,
    E_gm = E_gm*1e3,
    S_gm = S_gm*1e3,
    D_gm = D_gm*1e3,
    gs0 = gs0,
    a1 = a1,
    b1 = b1/1e3,
    gb = gb,
    Ncmin = Ncmin,
    angles = angles,
    sigma_soil = 2*sigma_soil, # Because it only reflects in one direction
    sky = sky,
    lat = toRadians(lat),
    Ca = Ca
)
end


# Variables will change during the season due to growth or acclimation
function variables(pars;
    # Parameters related to growth
    nleaf  = 8.0, # Leaf nitrogen content (g N/m2)
    wleaf = 200.0, # Leaf biomass as dry weight (g/m2)
    SLA = 300.0, # Specific leaf area (m2/g)
    # Variables related to distribution of nitrogen within canopy and leaf
    kN   = 0.4, # Extinction coefficient of leaf nitrogen
    f_Ncmn = 0.1, # Minimum fraction of N allocated to chlorophyll
    f_Ncmx = 0.3, # Maximum fraction of N allocated to chlorophyll
    pf_Nc  = 4.0, # Scaling constant for fraction of N allocated to chlorophyll
    f_Nrmn = 0.2, # Minimum fraction of N allocated to Rubisco
    f_Nrmx = 0.6, # Maximum fraction of N allocated to Rubisco
    pf_Nr  = 2.0 # Scaling constant for fraction of N allocated to Rubisco
)
    SLA = SLA/1e4
    # Constraint green LAI to avoid leaves with N < Nmin
    Ntotal = nleaf/14.0 # g/m2 -> mol/m2
    Nmin   = pars.Ncmin/SLA/14.0
    LAImax = (1/kN)*log(1 + kN*Ntotal/Nmin)
    GLAI = min(wleaf*SLA, LAImax)

    # Compute the leaf nitrogen content at top and bottom of the canopy
    Nlmax = kN*Ntotal/(1.0 - exp(-kN*GLAI))

(
    Nmin   = Nmin,
    GLAI = GLAI,
    Nlmax = Nlmax,
    kN = kN,
    f_Ncmn = f_Ncmn,
    f_Ncmx = f_Ncmx,
    pf_Nc  = pf_Nc,
    f_Nrmn = f_Nrmn,
    f_Nrmx = f_Nrmx,
    pf_Nr  = pf_Nr
)

end