# Dataframe and data types
using DataFrames
using DataFramesMeta
using CategoricalArrays
using Dates
using Statistics

# Data IO
using CSV
using JSONTables
using JSON
using GZip

# Testing
using BenchmarkTools
using GLMakie

function getColumnSettings(colSettingPath)
    return parseColumnSettings(JSON.parsefile(colSettingPath))
end

function parseColumnSettings(colSettings)
    parsed = deepcopy(colSettings) # So that we don't mutate the dict values
    for col in keys(parsed)
        parsed[col]["dtype"] = parseColDtype(colSettings[col]["dtype"])
        parsed[col]["procFun"] = parseColProcFun(colSettings[col]["proc"])
    end
    return parsed
end

function parseColDtype(dtypeString::String)
    # Parse the string specific of data type for a column
    if dtypeString == "str"
        return String31
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
        #dft1 = dateformat"y-m-dTH:M:S.s"
        dft2 = dateformat"yyyy-mm-ddTHH:MM:SS.ssssss"
        return col -> DateTime.(replace.(col, r"Z"=>""), dft2)
    else
        return nothing # We want uncaught cases to error out
    end
end
function parseColProcFun(nothing)
    # For the case where the json input is null, we don't apply any post-processing
    return identity
end


function firstNonMissingIfExists(arr)
    # Returns the first nonmissing value from an array if it exists, otherwise returns missing
    if length(collect(skipmissing(arr))) > 0
        return first(skipmissing(arr))
    else
        return missing
    end
end

function convertColumnTypes(df, colSettings)
    # by default the convert function doesn't allow missing values, so we create a version that does
    conversionFunWithMissing = passmissing(convert) 
    for colName in names(df)
        dtype = colSettings[colName]["dtype"]
        df[!,colName] = conversionFunWithMissing.(dtype, df[!,colName])
    end
    return df
end

function deduplicateUsingDataCols(df, colSettings)
    # Returns a dataframe that has exactly one row for each frame in the original data
    # In the case that the frames had missing data, this row may be a composite of multiple frames

    # Keep only the rows that have distinct values in the data columns (altitude, temperature, etc) ignoring packet metadata (uploader_callsign, uploaded_at, etc).
    df = unique(df,(Cols(col -> colSettings[String(col)]["kind"]=="data")))
    # split the dataframe into frames with only one row, and frames with multiple rows
    multipleEntryFrames= @chain df begin
        DataFramesMeta.groupby(:frame)
        subset(:frame=> f -> length(f) > 1, view=false, ungroup=false)
    end
    singleEntryFrames= @chain df begin
        DataFramesMeta.groupby(:frame)
        subset(:frame=> f -> length(f) == 1, view=false, ungroup=true)
    end
    # Construct the single most complete entry for each frame that has multiple distinct rows in the dataframe and add them to a list
    deduplicatedFrameRows = []
    for idx in eachindex(multipleEntryFrames)
        # take the subset of the original dataframe with rows for each particular frame
        framedf = multipleEntryFrames[idx]
        # Replace missing values with the first non-missing value if it exists
        transform!(framedf, DataFrames.All() .=> (x -> replace(x, missing => firstNonMissingIfExists(x))) => identity)
        # Overwrite the callsign to indicate that this is a composite frame
        transform!(framedf, Cols(:uploader_callsign) .=> (x -> String31("ZZZZZZ")) => identity)
        push!(deduplicatedFrameRows, framedf[1, DataFrames.All()])
    end
    df = vcat(singleEntryFrames, DataFrame(deduplicatedFrameRows))
    df = sort(df, :frame)
    return df
end

function applyPostProcessing(df, colSettings)
    for colName in names(df)
        f = colSettings[colName]["procFun"]
        df[!, colName] = f(df[!, colName])
    end
    return df
end

function dataframeFromJson(jsonPath, colSettings, convertTypes=true, deduplicate=true, postprocess=true)
    # Takes a json file from sondehub and returns a dataframe with columns selected and preprocessed according to colSettings 

    # Read in the gzipped json file
    df = DataFrame(jsontable(GZip.open(jsonPath)))
    # Keep relevant columns specified by colSettings
    df = select(df, (Cols(col -> colSettings[String(col)]["keep"]==true)))
    # Allow missing data in any column
    df = allowmissing(df)
    # Convert data types in case not correctly inferred
    if convertTypes
        df = convertColumnTypes(df, colSettings)
    end
    # Deduplicate rows intelligently, keeping one row for each frame while attempting to keep the most complete set of data possible for that frame
    if deduplicate
        df = deduplicateUsingDataCols(df, colSettings)
    end
    # Apply post-processing functions, if any, for each column
    if postprocess
        df = applyPostProcessing(df, colSettings)
    end
    return df
end

function csvFromDataframe(csvPath, df)
    # Save preprocessed dataframe to csv file
    CSV.write(csvPath, df, dateformat=dateformat"yyyy-mm-ddTHH:MM:SS.ssssss")
    return nothing
end

function csvFromJson(jsonPath, csvPath, colSettingPath)
    colSettings = getColumnSettings(colSettingPath)
    csvFromDataframe(csvPath, dataframeFromJson(jsonPath))
    return nothing
end

function dataframeFromCsv(csvPath::String, colSettingPath::String)
    colSettings = getColumnSettings(colSettingPath)
    df = dataframeFromCsv(csvPath, colSettings)
end
function dataframeFromCsv(csvPath::String, colSettings)
    # Load dataframe from csv file. This is assumed to be deduplicated
    # postprocessing functions must still be applied to convert from primitive json types to (for instance) DateTime
    # validate=false prevents CSV from throwing an error when there are columns in types that aren't found in the file
    df = CSV.read(csvPath, types=Dict(col => colSettings[col]["dtype"] for col in keys(colSettings)), DataFrame, validate=false)
    df = applyPostProcessing(df, colSettings)
    return df
end

jsonPath = "./sondehub/W1621037.json.gz"
csvPath = "./sondehub/W1621037_testSave.csv"
colSettingPath = joinpath([pwd(), "columnSettings", "hubColumnSettings.json"])
colSettings = getColumnSettings(colSettingPath)

df1 = dataframeFromJson(jsonPath, colSettings)
csvFromDataframe(csvPath, df1)
df2 = dataframeFromCsv(csvPath, colSettings)