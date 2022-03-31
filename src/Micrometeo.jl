####################################
############ Sky models ############
####################################

# Discretization of an uniform sky (UOC standard)
function uniform_sky(::Val{N}) where N
    Δθ = π/2/N
    θₗ  = SVector{N, Float64}(0.0:Δθ:π/2 - Δθ)
    θᵤ = SVector{N, Float64}(Δθ:Δθ:π/2)
    beta  = @. π/2 - (θₗ + θᵤ)/2
    f  = @. cos(θₗ)^2 - cos(θᵤ)^2
    return (beta = beta, f = f)
end


# Discretization of a standard overcast sky (SOC standard)
function standard_sky(::Val{N}) where N
    Δθ = π/2/N
    θₗ  = SVector{N, Float64}(0.0:Δθ:π/2 - Δθ)
    θᵤ = SVector{N, Float64}(Δθ:Δθ:π/2)
    beta  = @. π/2 - (θₗ + θᵤ)/2
    f  = @. ((4cos(θₗ) + 3)*cos(θₗ)^2 - (4cos(θᵤ) + 3)*cos(θᵤ)^2)/7
    return (beta = beta, f = f)
end

# Transforms an angle from hexadecimal degrees to radians
toRadians(angle) = angle*π/180.0
  
# Transforms an angle from radians into hexadecimnal degrees
toDegrees(angle) = angle*180.0/π
  
# Solar constant (J/m2/s)
const SCS = 1367.7
  
# Calculate the length of the diurnal part of a day based on latitude and day of year
function dayLength(lat, doy)::Float64
    # Declination angle of the Sun as seen from Earth.
    # The 10 is due to mismatch between winter solstice and the beginning of the Gregorian calendar
    dec = asin(sin(toRadians(-23.45))*cos(2π*(doy + 10.0)/365.0))
    # Calculate sunset angle with respect to solar noon, including corrections for extreme latitudes
    cosSunset = -tan(lat)*tan(dec)
    sunset = ifelse(cosSunset > 1.0, 0.0, ifelse(cosSunset < -1.0,  π, acos(cosSunset)))
    # dayLength in hours knowing that 1 hour = 15 degrees
    2.0*toDegrees(sunset)/15.0
end
  
# Returns time of day in solar time (hours)
function timeOfDay(time, DL)::Float64
    # Time of sunrise
    tsr = 12.0 - DL/2.0
    # diurnal time -> time of the day
    tod = tsr + time
    tod
end
  
# Calculate solar elevation angle
function calcBeta(lat, doy, time, DL)::Tuple{Float64, Float64}
    # Times of day at which we want to calculate solar radiation
    tod = timeOfDay(time, DL)
    # Angles of day in radians
    aod = toRadians((tod - 12.0)*15.0)
    # Declination angle of the sun
    dec = asin(sin(toRadians(-23.45))*cos(2π*(doy + 10.0)/365.0))
    # Solar elevation angle
    beta = asin(sin(lat)*sin(dec) + cos(lat)*cos(dec)*cos(aod))
    # Assume that when the sun is below the horizon there is no contribution to incoming solar radiation
    beta < 0.0 && (beta = 0.0)
    # Return the angle
    beta, dec
end

# Calculate solar azimuth angle
function calcAzimuth(beta, dec, lat, morning)::Float64
    Omega = acos(clamp((sin(dec) - sin(beta)*sin(lat))/(cos(beta)*cos(lat)), -1.0, 1.0))
    morning ? Omega : 2π - Omega
end


# Calculate solar radiation on the surface at a given time during the day
function solar_radiation(weather, time)::Float64
    # Calculate the solar angles
    beta, _ = calcBeta(weather.lat, weather.doy, time, weather.DL)
    # Solar irradiance on top of the atmosphere
    S0 = SCS*(1 + 0.033*cos(2π*(weather.doy - 1)/365.0))
    # Solar irradiance on Eath's surface (assume 75% transmissivity) in W/m2
    S0*weather.tau*sin(beta) 
end
  
# Calculate daily solar radiation outside of the atmosphere
function daily_solar_radiation(doy, lat, DL)::Float64
    # Extraterrestial solar constant (J/m2/s)
    S0 = SCS*(1 + 0.033*cos(2π*(doy - 1)/365.0))
    # Declination angle of the sun
    dec = asin(sin(toRadians(-23.45))*cos(2π*(doy + 10.0)/365.0))
    # Intermediate variables
    a = sin(lat)*sin(dec)
    b = cos(lat)*cos(dec)
    # Extraterrestial daily solar radiation (J/m2)
    S0*3600*(a*DL + 24/π*b*sqrt(1 - a^2/b^2))
end

  
# Calculate fraction of solar radiation that is diffuse based on empirical observations
function calculateDiffuse(Rsolar, weather, time)::Float64
    # Solar irradiance on top of the atmosphere project on horizontal plane
    S0 = SCS*(1 + 0.033*cos(2π*weather.doy/365.0))
    # Solar angles
    beta, _ = calcBeta(weather.lat, weather.doy, time, weather.DL)
    # Calculate Sg/S0
    SgS0 =  ifelse(beta == 0.0, 0.0, min(1.0, Rsolar/S0/sin(beta)))
    # Compute Sdf/Sg according to equation in Appendix at Spitters et al. (1986)
    R = 0.847 - 1.61*sin(beta) + 1.04*sin(beta)^2
    K = (1.47 - R)/1.66
    if SgS0 <= 0.22
        SdfSg = 1.0
    elseif SgS0 <= 0.35
        SdfSg = 1.0 - 6.4*(SgS0 - 0.22)^2
    elseif SgS0 <= K
        SdfSg = 1.47 - 1.66*SgS0
    else
        SdfSg = R
    end
    # Return the fractions
    SdfSg    
end
  

# Sine-exponential model by Ephrath et al. (1996, doi: 10.1016/0308-521X(95)00068-G)
# The flattening in the afternoon is not included
function daily_temperature(DL, Tmin, Tmax)
    # Precompute all variables (parameters treated as constants)
    LSH = 12.0 # becasue we are workin with solar time
    P = 1.0
    tsr = 12 - DL/2
    tss = 12 + DL/2
    NL = 24.0 - DL
    # Temperature at sunset
    Ts = Tmin + (Tmax - Tmin)*sin(π*(tss - LSH + DL/2)/(DL + 2P))
    tau = 4.0
    # Create anonymous function responsible for interpolation
    time -> begin 
        if time > tsr && time < tss
            Tmin + (Tmax - Tmin)*sin(π*(time - LSH + DL/2)/(DL + 2P))
        else
            time < tss && (time += 24.0)
            (Tmin - Ts*exp(-NL/tau) + (Ts  - Tmin)*exp(-(time - tss)/tau))/(1.0 - exp(-NL/tau))
        end

    end
end

# Calculate Ib0, Id0 and beta for a given day (clear day) at given time
function interpolate_meteo(weather, time)
    # Total solar radiation
    Rsolar = solar_radiation(weather, time)
    # Fraction that is diffuse
    fd = calculateDiffuse(Rsolar, weather, time)
    # Diffuse and direct (beam) radiation
    Ib0 = Rsolar*(1.0 - fd)
    Id0 = Rsolar*fd
    # Solar elevation angle
    beta, dec = calcBeta(weather.lat, weather.doy, time, weather.DL)
    # Solar azimuth angle
    ohm  = calcAzimuth(beta, dec, weather.lat, time < weather.DL/2)
    # Diurnal variation in temperature
    Tleaf = weather.Temp(time + weather.tsr)
    # VPD
    es = 0.61078*exp(17.269*Tleaf/(Tleaf + 237.3))
    VPD = max(0.0, es - weather.ea)
    # Return all relevant variables
    Environment(Ib0*0.45*4.57*1e-6, Id0*0.45*4.57*1e-6, beta, ohm,
                Tleaf + 273.15, VPD*1e3, weather.Ca*1e-6)
end

  
# Object that includes all daily weather variables as needed for interpolation
struct DailyWeather{T}
    lat::Float64
    doy::Int64
    tsr::Float64
    DL::Float64
    tau::Float64
    Temp::T
    ea::Float64
    Ca::Float64
end
function DailyWeather(lat::Float64, doy::Int64, tau::Float64, Tmin::Float64, 
                      Tmax::Float64, RH::Float64, Ca::Float64)
    DL = dayLength(lat, doy)
    tsr = 12.0 - DL/2
    Temp = daily_temperature(DL, Tmin, Tmax)
    Tsr = Temp(tsr)
    es = 0.61078*exp(17.269*Tsr/(Tsr + 237.3))
    ea = RH*es
    DailyWeather(lat, doy, tsr, DL, tau, Temp, ea, Ca)
end

# Object that includes all the relevant environmental variables for each time of the day
# Generated by interpolate_meteo
struct Environment
    Ib0::Float64
    Id0::Float64
    beta::Float64
    Omega::Float64
    Tleaf::Float64
    VPD::Float64
    Ca::Float64
end

