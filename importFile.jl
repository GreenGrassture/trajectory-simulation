using Profile
using PProf

if(!@isdefined typesImported)
    include("typeDefinitions.jl")
    typesImported = true
end

if (!(@isdefined setupComplete))||(setupComplete == false)
    include("physicalConstants.jl")
    include("trajectoryViewer.jl")
    include("balloonSimulator.jl")
    include("dataImporter.jl")
    setupComplete = true
end