module Elfel

# Elfel (C) 2020, Petr Krysl

include("utilities.jl")
include("LocalAssemblers.jl")
include("Assemblers.jl")
include("RefShapes.jl")
include("FElements.jl")
# include("FEFields.jl")
include("FESpaces.jl")
include("QPIterators.jl")
include("FEIterators.jl")

include("Exports.jl")

end # module
