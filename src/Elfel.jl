module Elfel

# Elfel (C) 2020, Petr Krysl

include("utilities.jl")
include("LocalAssemblers.jl")
include("Assemblers.jl")
include("RefShapes.jl")
include("FElements.jl")
include("FESpaces.jl")
include("QPIterators.jl")
include("FEIterators.jl")

# We can either use/import individual functions from Elfel like so:
# ```
# using Elfel.FESpaces: numberfreedofs!, numberdatadofs!, numberdofs!
# ```
# or we can bring into our context all exported symbols as
# ```
# using Elfel.Exports
# ```
include("Exports.jl")

end # module
