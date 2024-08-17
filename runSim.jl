using JSON
using Profile
using ProfileView
using FlameGraphs

include("importFile.jl")

# Contains physically "real" information like launch location, time, type of gas. Things that would apply to a real launch
launchParamPath = "launchParameters.json"
# Contains information pertaining to the simulation itself like maximum simulation time, initial position in simulation coordinates
simParamPath = "simParameters.json"


launchParams = getLaunchParameters(launchParamPath)
simParams = getSimParameters(simParamPath)

prob = constructSim(launchParams, simParams)
prob, sol, df = runSim(prob);
saveSim(prob, sol, df)

