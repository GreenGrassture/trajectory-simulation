# Until I expand the codebase into a proper module, this is where all the code loading happens.
# Scripts that actually *do* stuff should include this file. That way all the code loading and compilation isn't performed
# repeatedly.




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