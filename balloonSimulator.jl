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

using Dates
using TimeZones


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
    global R
    P = getAirPressure(u, p, t)
    T = getAirTemp(u, p, t)
    nGas = mGas/p[:liftGas].molarMass
    V = nGas*R*T/P
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

function dFun!(d, u, p, t)
    rx, ry, rz, vx, vy, vz, mGas, TGas = u
    ρt, Pt, Tt, μt = getAirParams(u, p, t) # air properties for this specific timestep
    uWind, vWind = getWindSpeed(u, p, t)
    fDrag = force_drag(u, p, t)
    fBuoyant = force_buoyancy(u, p, t)
    fGrav = force_gravity(p[:mBal], mGas, p[:mPay], g)
    ax, ay, az = (fDrag + fBuoyant + fGrav)/(mGas + p[:mBal] + p[:mPay])

    # dr/dt = v
    d.rx = vx
    d.ry = vy
    d.rz = vz
    # dv/dt = a
    d.vx = (uWind - vx)^3    # Should closely follow external wind velocity
    d.vy = (vWind - vy)^3
    d.vz = az
    # d/dt(mGas)
    d.mGas = 0.0 # contant for now
    # d/dt(TGas)
    d.TGas = (Tt - TGas)^3 # Should stay tightly locked onto external temperature
    return Nothing
end

function getGasParameters(gasJsonPath)
    # Retrieves the parameters of the gas (air or lift gas)
    gasJson = JSON.parsefile(gasJsonPath)
    return gasParams(gasJson["molar mass"], gasJson["contamination fraction"])
end

function getLaunchParameters(launchJsonPath)
    # Retrieves the launch parameters (latitude, longitude, time, balloon and payload mass...)
    launchJson = JSON.parsefile(launchJsonPath)
    return launchJson
end

function getSimParameters(simJsonPath)
    simJson = JSON.parsefile(simJsonPath)
    return simJson
end

function getLaunchLocation(launchParams)
    lat = launchParams["loc"]["lat"]
    lon = launchParams["loc"]["lon"]
    alt = launchParams["loc"]["alt"]
    loc = LLA(lat, lon, alt) 
    return loc
end

function getLaunchTime(launchParams)
    launchYear = launchParams["time"]["year"]
    launchMonth = launchParams["time"]["month"]
    launchDay = launchParams["time"]["day"]
    launchHour = launchParams["time"]["hour"]
    launchMinute = launchParams["time"]["minute"]
    launchSecond = launchParams["time"]["second"]
    launchTime = DateTime(launchYear, launchMonth, launchDay, launchHour, launchMinute, launchSecond)
    return launchTime
end

function constructSim(launchParams, simParams)
    launchPoint = getLaunchLocation(launchParams)
    launchTime = getLaunchTime(launchParams)

    # Define the transformation from local ENU (East, North, Up) coordinates in meters to global LLA (Lat, Lon, Alt) coordinates (altitude is in meters)
    loc2glob = LLAfromENU(launchPoint, grs80) # we use grs80 to be consistent with USGS data
    # Define the inverse transformation
    glob2loc = inv(loc2glob)

    # Initial position in local coordinates
    rx0, ry0, rz0 = simParams["pos"] # Height should be some nonzero value so the simulation doesn't immediately terminate
    # Initial velocity
    vx0, vy0, vz0 = simParams["vel"]

    ### Define lift gas parameters: initial mass, type of gas, etc.
    # Many of these values are sourced from Sóbester et al. 2014

    # ISA atmospheric data.  This function uses geometric altitude (ordinary meters)
    # Air parameters at launch. Units are kg/m^3, Pa, K, ???
    ρ0, P0, T0, μ0 = ISAdata(launchPoint.alt)

    # Parameters of the lift gas. In this case hydrogen
    gasPure = getGasParameters("gas_hydrogen.json")
    # Parameters of ambient air
    air = getGasParameters("gas_air.json")

    # Weight the molar mass of our lift gas toward that of ambient air to account for impurities
    molarMassGas = (gasPure.molarMass)*(1-gasPure.contamFrac) + gasPure.contamFrac*air.molarMass

    gasLift = gasParams(molarMassGas, 0.0)

    densityGasInitial = P0*molarMassGas/(R*T0)

    radiusInitial = 0.9 # m - initial balloon radius. 
    volumeInitial = (4/3)*pi*(radiusInitial^3) # m^3 - initial balloon volume
    # the mass of lifting gas is assumed to be known.  In this case case from a known (approximate) initial radius
    mGas0 = volumeInitial*densityGasInitial  # kg mass of the lifting gas
    mBal =  1.5 # kg - mass of the balloon
    mPay = 0.8 # kg - mass of the payload

    # Set up the callback that causes the balloon to burst when it becomes too large
    burstDiameter = launchParams["burst diameter"] # m

    u0 = ComponentArray(rx=rx0, ry=ry0, rz=rz0, vx=vx0, vy=vy0, vz=vz0, mGas=mGas0, TGas=T0)

    # Initial and final times in seconds
    t0::Float64 = simParams["initial time"]
    tFinal::Float64 = simParams["max sim duration"]
    tSpan = [t0, tFinal]

    p = Dict(
    :isBurst=>false, 
    :mBal=>mBal, 
    :mPay=>mPay, 
    :mGas=>mGas0,
    :liftGas=>gasLift, 
    :burstDiameter=>burstDiameter, 
    :l2g=>loc2glob,
    :g2l=>glob2loc,
    :launchLoc=>launchPoint,
    :launchTime=>launchTime)
    # Turn the dict into a namedtuple so it is immutable
    p = NamedTuple(pairs(p))

    cb_burst = ContinuousCallback(cb_burst_condition, cb_burst_affect!)
    cb_negativeAltitude = ContinuousCallback(cb_negativeAltitude_condition, nothing, cb_negativeAltitude_affect!)
    cb_all = CallbackSet(cb_burst, cb_negativeAltitude)

    # Attach symbolic names 
    f! = ODEFunction(dFun!, syms=[Symbol(s) for s in ComponentArrays.labels(u0)])
    prob = ODEProblem(f!, u0, tSpan, p, callback=cb_all)
    return prob
end

function runSim(prob)
    sol = solve(prob)
    # Since we atteched the symbolic names, we can just transform this into a dataframe with the appropriate column names like this:
    df = DataFrame(sol)
    # Convert to DateTim
    return (prob, sol, df)
end

function saveSim(prob, sol, df)
    df = rename(df, ("timestamp"=>"datetime"))
    # The "native" unit of the ODE solution is seconds. The choice to convert to nanoseconds is somewhat arbitrary
    # but in theory it should preserve the most significant figures when rounding to an integer for type conversion
    df[!, :datetime] = prob.p[:launchTime] + Nanosecond.(round.(df[!, :datetime]*1e9))

    # Trajectory in global LLA coordinates
    llas = Base.stack((lla -> [lla.lat, lla.lon, lla.alt]).([prob.p[:l2g](ENU(Array(row[[:rx, :ry, :rz]]))) for row in eachrow(df)]))'
    df[!,[:rx, :ry, :rz]] = llas
    rename!(df, [:rx, :ry, :rz] .=> [:lat, :lon, :alt])

    simFolder = "./simulations/"
    simFile = "sol.csv"
    CSV.write(joinpath(simFolder,simFile), df, dateformat=dateformat"yyyy-mm-ddTHH:MM:SS.ssssss")
    #idxRz = label2index(u0, "rz")
end