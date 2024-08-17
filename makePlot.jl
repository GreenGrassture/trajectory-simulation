include("importFile.jl")

# Import our test trajectories
simFolder = "./simulations/"
simFile = "sol.csv"
simPath = joinpath(simFolder, simFile)
dfSim = dataframeFromCsv(simPath)

hubFolder = "./sondehub/"
hubFile = "W1621037_testSave.csv" #"V4630041.csv"
hubPath = joinpath(hubFolder, hubFile)
colSettingPath = joinpath(cd(pwd, ".."), "columnSettings.json")
dfHub = 

geoDir = "geography/one arc sec/"
geoBig = importGeography(geoDir)

fig = plotTrajectories([dfSim, dfHub])
display(fig)