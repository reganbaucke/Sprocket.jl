using JuMP

mutable struct Problem
    vars::Set{Variable}
    model::JuMP.Model
    c_oracle::Function
    m_oracle::Function
    domain
end

function Problem()
    variables = Set{Variable}()
    model = JuMP.Model()
    c_oracle(x...) = nothing
    m_oracle(x...) = nothing
    return Problem(variables,model,c_oracle,m_oracle,())
end


function Sprocket.add_variable(p::Problem,var::Variable)
    #check that the variable name is not already in the model
    @assert !in(var.name,map(x->x.name,collect(p.vars)))
    #
    push!(p.vars,var)
end

function set_domain(p::Problem,domain)
    #verify that the point has every random variable in it and nothing more
    for point in domain
        @assert filter(x -> x.type == Random(),p.vars) == get_vars(domain[1]) == get_vars(domain[2])
    end
    p.domain = domain
    return nothing
end
