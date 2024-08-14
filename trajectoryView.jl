using Statistics

using CSV
using DataFrames
using Dates

# Handling geotiff data
using CoordinateTransformations
using Rotations
using StaticArrays
using Geodesy
using TiffImages
using Rasters
using ArchGDAL

# 3d visualization
using GLMakie
using Mux

#using Images
#using FileIO

# Defining the origin of the coordinate system for the trajectory
lat = 37.7195
lon = -122.4191
alt = 81
launchPoint = LLA(lat, lon, alt) # Field of Dogs, McLaren Park, SF
# Define the transformation from local East, North, Up coordinates (in meters) to global Lat, Lon, Alt coordinates
loc2glob = LLAfromENU(launchPoint, grs80) # we use grs80 because that's what the USGS uses
glob2loc = inv(loc2glob)


function importSimData(simFolder, simFile)
    df = DataFrame(CSV.File(joinpath(simFolder,simFile)))

    rx = df[!, :rx]
    ry = df[!, :ry]
    rz = df[!, :rz]
    t = df[!, :t]
    # Trajectory in global coordinates
    rsGlobal = [loc2glob(ENU(rx[i], ry[i], rz[i])) for i in eachindex(rx)]
    rLats = [r.lat for r in rsGlobal]
    rLons = [r.lon for r in rsGlobal]
    rAlts = [r.alt for r in rsGlobal]
    return (t, rLats, rLons, rAlts)
end

function importHubData(hubFolder, hubFile)
    df = DataFrame(CSV.File(hubFolder*hubFile))
    rLats = df[!,:lat]
    rLons = df[!,:lon]
    rAlts = df[!,:alt]
    datetimes = DateTime.(String.(df[!,:datetime]), dateformat"yyyy-mm-ddTHH:MM:SS.ssssssZ")
    # extract the milliseconds since launch, then convert to seconds. We can't get the seconds
    # directly because of some issue with rounding when converting
    t = (Dates.value.(datetimes .- datetimes[1]))./1000.0
    
    return (t, rLats, rLons, rAlts)
end

# Import our test trajectory
simFolder = "./simulations/"
simFile = "sol.csv"
simTime, simLats, simLons, simAlts = importSimData(simFolder, simFile)

hubFolder = "./sondehub/"
hubFile = "W1621037.csv" #"V4630041.csv"
hubTime, hubLats, hubLons, hubAlts = importHubData(hubFolder, hubFile)

geoDir = "geography/one arc sec/"
tifPaths = readdir(geoDir)

rasters = [Raster(geoDir*path) for path in tifPaths]
rasters = [resample(rast, size=(500, 500)) for rast in rasters]

geoBig = Rasters.mosaic(mean, rasters, atol=0.001)

#=
tifPathN38W123 = "geography/one arc sec/USGS_1_n38w123_20240313.tif"
rasN38W123 = Raster(tifPathN38W123, lazy=false)
rasN38W123 = resample(rasN38W123, size=(1000, 1000), method=:average)
=#

GLMakie.activate!
set_theme!(backgroundcolor = :lightskyblue2)
aspect=(1, 1, 1)
perspectiveness=0.2
fig= Figure(; size=(800, 800))
ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)

xlims!(ax1, (-124, -119))
ylims!(ax1, (36, 39))
zlims!(ax1, (-10.0, 80000.0))
#contourf!(ax, geoSmall, colormap=:batlow)
surface!(ax1, geoBig, colormap=:oleron, nan_color=:red, colorrange=(-1499,1500))

#surface!(ax1, rasN38W122, colormap=:winter, nan_color=:red, colorrange=(-1,1500))
#surface!(ax1, rasN38W123, colormap=:winter, nan_color=:red, colorrange=(-1,1500))

lines!(ax1, simLons, simLats, simAlts, color=simTime, linewidth=2, colormap=:winter)
lines!(ax1, hubLons, hubLats, hubAlts, color=hubTime, linewidth=2, colormap=:lajolla)

#cam = Makie.Camera3D(ax1.scene, projectiontype = Makie.Perspective, reset = Keyboard.left_control, center=false, reposition_button=Keyboard.left_alt)
display(fig)