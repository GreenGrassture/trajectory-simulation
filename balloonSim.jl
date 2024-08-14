# An attempt at modeling the volume of a balloon and related variables (pressure, buoyant force etc)

using LinearAlgebra
using Statistics
using DifferentialEquations
#using ModelingToolkit
using ISAData
#using Plots
using Geodesy

using Parameters
using ComponentArrays

using DataFrames
using CSV

using BenchmarkTools




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
    # returns ρ, P, T, μ
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    return ISAdata(rz) 
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
    rx, ry, rz, vx, vy, vz, mGas, TGas = u

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
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    isBurst = p[:isBurst]
    if isBurst
        area = 2.0
    else
        radius = getBalloonRadius(getBalloonVolume(u, p, t))
        area = pi*radius^2
    end
    return area
end

function force_buoyancy(u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    isBurst = p[:isBurst]

    if isBurst
        fBuoyant = [0.0; 0.0; 0.0]
    else
        V = getBalloonVolume(u, p, t)
        ρAir = getAirDensity(u, p, t)
        fBuoyant = [0.0; 0.0; V*ρAir*g]
    end
    return fBuoyant
end

function force_drag(u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    ρAir = getAirDensity(u, p, t)
    A = getProjectedArea(u, p, t)
    Re = getReynoldsNumber(u, p, t)
    # If the balloon has burst, the drag coefficient is determined by the parachute (mostly)
    if p[:isBurst]
        Cd = 0.7
    else
        Cd = getBalloonDragCoeff(Re)
    end

    fDrag = [0.0; 0.0; -(1/2)*ρAir*A*Cd*(vz^2)*sign(vz)]
    return fDrag
end

function force_gravity(mBal, mGas, mPay, g)
    return [0.0; 0.0; -g*(mBal + mGas + mPay)]
end

function getWindSpeed(u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    h = 10000.0
    f = 3600.0
    return  [sin(rz/h), cos(rz/h)] + 10.1*[sin(t/f), sin(t/f)] #[10.0, 0.0]
end

function getReynoldsNumber(L, v, μ, ρ)
    return ρ*v*L/μ
end

function getReynoldsNumber(u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    V = getBalloonVolume(u, p, t)
    D = 2*getBalloonRadius(V)
    ρ, P, T, μ = getAirParams(u, p, t)
    return getReynoldsNumber(D, vz, μ, ρ)
end

function getBalloonDragCoeff(Re)
    CdLaminar = 0.425 
    CdTurbulent = 0.225
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
ρ0, P0, T0, μ0 = ISAdata(0.0)

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

# Set up the callback that causes the balloon to burst when it becomes too large
burstDiameter = 9.45 # m
function cb_burst_condition(u, t, integrator)
    isBurst = integrator.p[:isBurst]
    # Checking to see if callback has returned true before. If so, this prevents it from running again
    # https://discourse.julialang.org/t/differentialequations-callback-which-is-only-triggered-the-first-time-condition-is-satisfied/70216/3
    if isBurst
        return 2.0
    else
        burstDiameter = integrator.p[:burstDiameter]
        rx, ry, rz, vx, vy, vz, mGas, TGas = u
        V = getBalloonVolume(u, integrator.p, t)
        R = getBalloonRadius(V)
        balloonDiameter = 2*R
        return burstDiameter - balloonDiameter
    end
end
function cb_burst_affect!(integrator)
    # p is a namedtuple, which is immutable, so we convert to a dict, change one key,
    # then convert back.  This is probably fine since it should only happen once during a run
    p = Dict(pairs(integrator.p))
    p[:isBurst] = true
    p = NamedTuple(pairs(p))
    integrator.p = p
    return Nothing
end

# Callback to terminate the simulation at zero altitude
function cb_negativeAltitude_condition(u, t, integrator)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    return rz
end 

function cb_negativeAltitude_affect!(integrator)
    return terminate!(integrator)
end

cb_burst = ContinuousCallback(cb_burst_condition, cb_burst_affect!)
cb_negativeAltitude = ContinuousCallback(cb_negativeAltitude_condition, nothing, cb_negativeAltitude_affect!)
cb_all = CallbackSet(cb_burst, cb_negativeAltitude)

# 37.762380, -122.418405 initial location NB
# Initial position
rx0 = 0.0
ry0 = 0.0
rz0 = 1.0 # Height some nonzero value so the simulation doesn't immediately terminate

# Initial velocity
vx0 = 0.0
vy0 = 0.0
vz0 = 0.0

mGas0 = mGas
TGas0 = temperatureInitialAmbient

u0 = ComponentArray(rx=rx0, ry=ry0, rz=rz0, vx=vx0, vy=vy0, vz=vz0, mGas=mGas0, TGas=TGas0)

# Initial and final times in seconds
t0 = 0.0
tFinal = 10000.0
tSpan = [0.0, 10000.0]

p = Dict(:isBurst=>false, :mBal=>mBal, :mPay=>mPay, :mGas=>mGas, :burstDiameter=>burstDiameter)
p = NamedTuple([pair for pair in p])

function dFun!(d, u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    ρt, Pt, Tt, μt = getAirParams(u, p, t) # air parameters for this specific timestep
    uWind, vWind = getWindSpeed(u, p, t)
    fDrag = force_drag(u, p, t)
    fBuoyant = force_buoyancy(u, p, t)
    fGrav = force_gravity(mBal, mGas, mPay, g)
    ax, ay, az = (fDrag + fBuoyant + fGrav)/(mGas + mBal + mPay)

    # dr/dt = v
    d.rx = vx
    d.ry = vy
    d.rz = vz

    # dv/dt = a
    # Should closely follow external wind velocity
    d.vx = (uWind - vx)^3
    d.vy = (vWind - vy)^3
    d.vz = az

    # d/dt mGas
    d.mGas = 0.0 # contant for now
    # d/dt TGas
    d.TGas = (Tt - TGas)^3 # Should stay tightly locked onto external temperature
    return Nothing
end
# This will eventually be used when we have to upgrade to a proper DAE but for now we can fake it 
#M = Diagonal([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0])

# Attach symbolic names 
f! = ODEFunction(dFun!, syms=[Symbol(s) for s in labels(u0)])
prob = ODEProblem(f!, u0, tSpan, p, callback=cb_all)
sol = solve(prob)

df = DataFrame(sol)
df = rename(df, ("timestamp"=>"t"))
simFolder = "./simulations/"
simFile = "sol.csv"
CSV.write(simFolder*simFile, df)
#idxRz = label2index(u0, "rz")


