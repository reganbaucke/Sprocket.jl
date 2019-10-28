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


# define map for option values, how is this not in the standard library?
function Base.map(f, thing::Some{T}) where T
  thing.value |> f |> Some
end

function Base.map(f, thing::Nothing)
  nothing
end

function flat_map(f, thing::Some{T}) where T
  thing |> get |> f
end

function flat_map(f, thing::Nothing) where T
  nothing
end

function >>(thing,f)
  flat_map(f,thing)
end
