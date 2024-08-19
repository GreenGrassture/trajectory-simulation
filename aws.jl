# Dataframe and data types
using DataFrames
using DataFramesMeta
using CategoricalArrays
using Dates
using Statistics

# Data IO
using FilePathsBase
using CSV
using JSONTables
using JSON
using GZip

# AWS
using AWS
using AWSCore
using AWSS3

# Testing
using BenchmarkTools

using HTTP
using FilePathsBase
using FilePathsBase: /

using StringEncodings


function getTimeseriesFromSerial(serialNum, folderPath)
    hubLoc = "http://sondehub-history.s3.amazonaws.com/serial/"
    extn = ".json.gz"
    loc = reduce(*, [hubLoc, serialNum, extn])
    response = HTTP.get(loc)
    meta = response.headers
    decoded = GZip.read(decode(response.body, "UTF-8"))
    #JSON.write(joinpath([folderPath, serialNum, extn]), decoded)
    return decoded
end

folderPath = joinpath(pwd(), "sondehub")
serial = "W1621035"

dec = getTimeseriesFromSerial(serial, folderPath)
