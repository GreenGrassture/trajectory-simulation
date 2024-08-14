# Dataframe and data types
using DataFrames
using DataFramesMeta
using CategoricalArrays
using Dates

# Data IO
using CSV
using JSONTables
using JSON
using GZip

# Testing
using BenchmarkTools
using GLMakie

function parseColumnSettings(colSettings)
    parsed = deepcopy(colSettings) # So that we don't mutate the dict values
    for col in keys(parsed)
        parsed[col]["dtype"] = parseColDtype(colSettings[col]["dtype"])
        parsed[col]["procFun"] = parseColProcFun(colSettings[col]["dtype"])
    end
    return parsed
end

function parseColDtype(dtypeString::String)
    # Parse the string specific of data type for a column
    if dtypeString == "str"
        return String
    elseif dtypeString == "float64"
        return Float64
    elseif dtypeString == "int64"
        return Int64
    end
end

function parseColProcFun(procFunString::String)
    # Parse the string specification of postprocessing for a column
    if procFunString == "categorical"
        return col -> categorical(col, compress=true)
    elseif procFunString == "datetime"
        return col -> DateTime.(col, dateformat"yyyy-mm-ddTHH:MM:SS.ssssssZ")
    else
        return identity
    end
end

function colDtype(col, colSettings)
    return
end

function colKeep(col, colSettings)
    return
end

df = DataFrame(jsontable(GZip.open("./sondehub/W1621037.json.gz")))
df = allowmissing(df)
colSettingPath = joinpath(cd(pwd, ".."), "columnSettings.json")
colSettings = parseColumnSettings(JSON.parsefile(colSettingPath))

# Remove ignored columns specified by colSettings
dropCols = []
for col in names(df)
    if colSettings[col]["keep"] == false
        push!(dropCols, Symbol(col))
    end
end
df = df[:, Not(Cols(dropCols))]

# by default the convert function doesn't allow missing values
conversionFunWithMissing = passmissing(convert) 
for colName in names(df)
    dtype = colSettings[colName]["dtype"]
    df[!,colName] = conversionFunWithMissing.(dtype, df[!,colName])
end

# Keep only the rows that have distinct values in the data columns (altitude, temperature, etc) ignoring packet metadata like uploader etc.
dataCols = []
for col in names(df)
    if colSettings[col]["kind"] == "data"
        push!(dataCols, Symbol(col))
    end
end
df = unique(df,(Cols(dataCols)))

gd = groupby(df, Cols(:frame))

mult = @chain df begin
        groupby(:frame)
        filter(:frame=> f -> length(f) > 1, _)
end