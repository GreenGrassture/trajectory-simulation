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

function convertDatetimesToTimeSinceLaunch(datetimes)
    # Date.value is in milliseconds
    timeSinceLaunch = (Dates.value.(datetimes .- datetimes[1]))./1000.0
    return timeSinceLaunch
end

function importGeography(geoDir)
    tifPaths = readdir(geoDir)
    rasters = [Raster(geoDir*path) for path in tifPaths]
    rasters = [resample(rast, size=(500, 500)) for rast in rasters]
    geoBig = Rasters.mosaic(mean, rasters, atol=0.001)
    return geoBig
end

function plotTrajectories(trajectoryDfs, geoBig)
    GLMakie.activate!
    set_theme!(backgroundcolor = :lightskyblue2)
    aspect=(1, 1, 1)
    perspectiveness=0.2
    fig= Figure(; size=(800, 800))
    ax1 = Axis3(fig[1, 1]; aspect, perspectiveness)

    xlims!(ax1, (-124, -119))
    ylims!(ax1, (36, 39))
    zlims!(ax1, (-10.0, 80000.0))
    surface!(ax1, geoBig, colormap=:oleron, nan_color=:red, colorrange=(-1499,1500))
    for df in trajectoryDfs
        timeSeconds = convertDatetimesToTimeSinceLaunch(df[!,:datetime])
        lines!(ax1, df[!,:lon], df[!,:lat], df[!,:alt], 
               color=timeSeconds, linewidth=2, colormap=:lajolla)
    end
    #cam = Makie.Camera3D(ax1.scene, projectiontype = Makie.Perspective, reset = Keyboard.left_control, center=false, reposition_button=Keyboard.left_alt)
    return fig
end
