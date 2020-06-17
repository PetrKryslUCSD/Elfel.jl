module AggregateFEIterators

using StaticArrays
using ..FEIterators: FEIterator
using ..Utilities: @_check

#= TODO is it more natural to have access to the geometry from the font element space or from the iterator? =#
struct AggregateFEIterator
    iterators::Vector{FEIterator}
    nel::Int64

    function AggregateFEIterator(iterators)
        nel = Base.length(iterators[1])
        for i in iterators
            @_check Base.length(i)  == nel
        end
        return new(iterators, nel)
    end
end



function Base.iterate(it::AggregateFEIterator, state = 1)
    if state > it.nel
        return nothing
    else
        return (_update!(it, state), state+1)
    end
end
Base.length(it::AggregateFEIterator)  = it.nel

function _update!(it::AggregateFEIterator, state)
    for i in it.iterators
        Base.iterate(i, state)
    end
    return it
end

iterator(it::AggregateFEIterator, j) = it.iterators[j]
    

end
