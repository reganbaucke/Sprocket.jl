using JuMP

mutable struct Problem
    vars::Set{Variable}
    model::JuMP.Model
    c_oracle
    m_oracle
end

function Problem()
    variables = Set{Variable}()
    model = JuMP.Model()
    c_oracle(x...) = nothing
    m_oracle(x...) = nothing
    return Problem(variables,model,c_oracle,m_oracle)
end


function Sprocket.add_variable(p::Problem,var::Variable)
    @assert !in(var.name,map(x->x.name,collect(p.vars)))
    push!(p.vars,var)
end
