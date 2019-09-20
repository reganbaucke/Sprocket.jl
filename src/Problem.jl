mutable struct Problem
    vars
    model
    c_oracle
    m_oracle
end

function Problem()
    variables = Dict()
    model = JuMP.Model()
    c_oracle(x...) = nothing
    m_oracle(x...) = nothing
    return Problem(variables,model,c_oracle,m_oracle)
end
