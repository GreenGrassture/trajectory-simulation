# An attempt at modeling the volume of a balloon and related variables (pressure, buoyant force etc)

using DifferentialEquations
using ISAData
using Plots



# https://en.wikipedia.org/wiki/International_Standard_Atmosphere
# See notes under table on above page for explanation of calculation
function heightToGeometric(h)
    # h: geopotential height in geopotential meters
    r = 6356766
    z = r*h/(r - h)
    return z
end

function heightToGeopotential(z)
    # z: geometric height in meters
    r = 6356766
    h = r*z/(r + z)
    return h
end

# These functions will likely get more complicated when I integrate information from weather data, but for now they just use ISA data
function getAirParams(u, p, t)
    z, _ = u
    return ISAdata(z)
end

function getAirDensity(u, p, t)
    return getAirParams(u, p, t)[1]
end

function getAirPressure(u, p, t)
    return getAirParams(u, p, t)[2]
end

function getAirTemp(u, p, t)
    return getAirParams(u, p, t)[3]
end

function getGasViscosity(u, p, t)
    return getAirParams(u, p, t)[4]
end

function getGasDensity(u, p, t)
    # Computes the density of the balloon lift gas
    V = getBalloonVolume(u, p, t)
    return mGas/V
end

function getBalloonVolume(u, p, t)
    # For now, we're assuming a lot about the balloon:
    # 1) The internal temperature is uniform, and the same as the external temperature
    # 2) The pressure is the same as the external pressure
    # 3) The lift gas follows the ideal gas law
    z, vz = u

    P = getAirPressure(u, p, t)
    T = getAirTemp(u, p, t)
    V = mGas*RSpec*T/P
    return V
end

function getBalloonRadius(vBalloon)
    # Assuming a spherical balloon, get the radius from the volume
    return (vBalloon*(3/(4*pi)))^(1/3)
end

function getProjectedArea(u, p, t)
    z, vz = u
    isBurst = p["isBurst"]
    mGas = p["mGas"]
    if isBurst
        area = 2.0
    else
        radius = getBalloonRadius(getBalloonVolume(u, p, t))
        area = pi*radius^2
    end
    return area
end

function force_buoyancy(u, p, t)
    z, vz = u
    mGas = p["mGas"]
    isBurst = p["isBurst"]

    if isBurst
        fBuoyant = 0.0
    else
        V = getBalloonVolume(u, p, t)
        ρAir = getAirDensity(u, p, t)
        fBuoyant = V*ρAir*g
    end
    return fBuoyant
end

function force_drag(u, p, t)
    z, vz = u
    ρAir = getAirDensity(u, p, t)
    A = getProjectedArea(u, p, t)
    Re = getReynoldsNumber(u, p, t)
    # If the balloon has burst, the drag coefficient is determined by the parachute (mostly)
    if p["isBurst"]
        Cd = 0.8
    else
        Cd = getDragCoeff(Re)
    end

    fDrag = -(1/2)*ρAir*A*Cd*(vz^2)*sign(vz)
    return fDrag
end

function force_gravity(mBal, mGas, mPay, g)
    return -g*(mBal + mGas + mPay)
end

function getReynoldsNumber(L, v, μ, ρ)
    return ρ*v*L/μ
end

function getReynoldsNumber(u, p, t)
    z, vz = u
    V = getBalloonVolume(u, p, t)
    D = 2*getBalloonRadius(V)
    ρ, P, T, μ = getAirParams(u, p, t)
    return getReynoldsNumber(D, vz, μ, ρ)
end

function getDragCoeff(Re)
    CdLaminar = 0.225
    CdTurbulent = 0.425
    Re1 = 3.296e5
    ΔRe = 0.363e5
    Re2 = Re1 + ΔRe

    if Re < Re1
        Cd = CdLaminar
    elseif Re < Re2
        Cd = CdLaminar - (CdLaminar - CdTurbulent)*(Re - Re1)/ΔRe
    else
        Cd = CdTurbulent
    end
    return Cd
end

# Flag for whether we are using hydrogen or helium
usingH2 = true

# ISA atmospheric data at sea level.  This function uses geometric altitude (ordinary meters)
ρ0, P0, T0, μ0 = ISAdata(0)

pressureInitialAmbient = P0  # Pa - Air pressure at launch. For now, assume one standard atmosphere of pressure
temperatureInitialAmbient = T0 # K - Air temperature at launch 15°C

molarMassHydrogen =  0.00201 # kg∕mol - molar mass of hydrogen gas
molarMassHelium = 0.004 # kg∕mol

fracAirHelium = 0.055 # fraction of lift gas contamination with ambient air for helium
fracAirHydrogen = 0.015 # fraction of lift gas contamination with ambient air for hydrogen
fracAir = usingH2 ? fracAirHydrogen : fracAirHelium # as before, we are using hydrogen in our case

molarMassAir = 0.02896 #kg∕mol - molar mass of ambient air
# Weight the molar mass of our lift gas toward that of ambient air to account for impurities
molarMassGas = (usingH2 ? molarMassHydrogen : molarMassHelium)*(1-fracAir) + fracAir*molarMassAir

R = 8.31447 # m^3*Pa/(mol*K) - ideal gas constant
RSpec = R/molarMassGas # Specific ideal gas constant for our lift gas for equation of the form PV = mRT

densityAirInitial = pressureInitialAmbient*molarMassAir/(R*temperatureInitialAmbient)
densityGasInitial = pressureInitialAmbient*molarMassGas/(R*temperatureInitialAmbient)

radiusInitial = 0.9 # m - initial balloon radius. 
volumeInitial = (4/3)*pi*(radiusInitial^3) # m^3 - initial balloon volume
# the mass of lifting gas is assumed to be known.  In this case case from a known (approximate) initial radius
mGas = volumeInitial*densityGasInitial  # kg mass of the lifting gas
mBal =  1.5 # kg - mass of the balloon
mPay = 0.8 # kg - mass of the payload

g = 9.81 # m/s^2

function dFun(u, p, t)
    z, vz = u
    fDrag = force_drag(u, p, t)
    fBuoyant = force_buoyancy(u, p, t)
    fGrav = force_gravity(mBal, mGas, mPay, g)
    az = (fDrag + fBuoyant + fGrav)/(mGas + mBal + mPay)
    return [vz, az]
end

# Set up the callback that causes the balloon to burst when it becomes too large
burstDiameter = 9.45 # m
function cb_burst_condition(u, t, integrator)
    isBurst = integrator.p["isBurst"]
    # Checking to see if callback has returned true before. If so, this prevents it from running again
    # https://discourse.julialang.org/t/differentialequations-callback-which-is-only-triggered-the-first-time-condition-is-satisfied/70216/3
    if isBurst
        return 2.0
    else
        burstDiameter = integrator.p["burstDiameter"]
        z, vZ = u
        V = getBalloonVolume(u, integrator.p, t)
        R = getBalloonRadius(V)
        balloonDiameter = 2*R
        return burstDiameter - balloonDiameter
    end
end
function cb_burst_affect!(integrator)
    integrator.p["isBurst"] = true
end

# Callback to terminate the simulation at zero altitude
function cb_negativeAltitude_condition(u, t, integrator)
    z, vz = u
    return z
end 

function cb_negativeAltitude_affect!(integrator)
    return terminate!(integrator)
end

cb_burst = ContinuousCallback(cb_burst_condition, cb_burst_affect!)
cb_negativeAltitude = ContinuousCallback(cb_negativeAltitude_condition, nothing, cb_negativeAltitude_affect!)
cb_all = CallbackSet(cb_burst, cb_negativeAltitude)
zInitial = 1.0
zDotInitial = 0.0
u0 = [zInitial, zDotInitial]
tSpan = [0.0, 10000.0]
isBurst = false

p = Dict("isBurst"=>isBurst, "mBal"=>mBal, "mPay"=>mPay, "mGas"=>mGas, "burstDiameter"=>burstDiameter)


prob = ODEProblem(dFun, u0, tSpan, p, callback=cb_all)
sol = solve(prob)

plot(sol)