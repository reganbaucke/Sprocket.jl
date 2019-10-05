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
