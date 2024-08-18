include("importFile.jl")

colSettingFolder = joinpath([pwd(), "columnSettings"])
hubColSettingPath = joinpath([colSettingFolder, "hubColumnSettings.json"])
simColSettingPath = joinpath([colSettingFolder, "simColumnSettings.json"])

# Import our test trajectories
simFolder = joinpath([pwd(),"simulations"])
simFile = "sol.csv"
simPath = joinpath([simFolder, simFile])
dfSim = dataframeFromCsv(simPath, simColSettingPath)

hubFolder = joinpath([pwd(), "./sondehub/"])
hubFile = "W1621037_testSave.csv" #"V4630041.csv"
hubPath = joinpath(hubFolder, hubFile)
dfHub = dataframeFromCsv(hubPath, hubColSettingPath)

geoDir = "geography/one arc sec/"
geoBig = importGeography(geoDir)

fig = plotTrajectories([dfSim, dfHub], geoBig)
display(fig)