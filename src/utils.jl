############
# Utility functions and miscillaneous helpers
############
using JuMP

#### delete whole array of variables from a jump model
function JuMP.delete(m::JuMP.Model,vars::AbstractArray)
    for var in vars
        JuMP.delete(m,var)
    end
end

function JuMP.size(::JuMP.VariableRef)
    return ()
end


# Option values for error handling in a functional way
struct Some{T}
    value::T
end

const Option{T} = Union{Some{T},Nothing}

function get(thing:: Some{T}) where T
    thing.value
end

function Base.map(f,thing::Some{T}) where T
    thing |> get |> f |> Some
end

function Base.map(f,thing::Nothing)
    nothing
end

function flat_map(f,thing::Some{T}) where T
    thing |> get |> f
end

function flat_map(f,thing::Nothing) where T
    nothing
end

function some_function(x::Number)
    if x <= 5
        Some(x+6)
    else
        nothing
    end
end

function >>(thing::Option,f)
    flat_map(f,thing)
end
