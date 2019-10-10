module Sprocket
# include("./Point.jl")
include("./Variable.jl")
include("./Cut.jl")
include("./Problem.jl")
include("./Bound.jl")
include("./Criteria.jl")
include("./utils.jl")

# the main solving routine

# function Sprocket.solve(algorithm,prob::Problem,my_crit::Criteria)
# 	(init,iterate,check) = algorithm
# 	states = []
# 	push!(states,init(prob))
# 	while !check(states[end])
# 		push!(states,iter(states[end]))
# 	end
# 	return states
# end

function Sprocket.solve(algorithm,prob,crit)
	(init,iterate,check) = algorithm
	states = []
	push!(states,init(prob))
	# while !check(crit,states[end])
	for i = 1:30
		push!(states,iterate(states[end],prob))
	end
	# end
	return states
end

# ####
# # this is the interface that I want to present to the user
# ####
# function build_problem()
# 	# problem = Sprocket.Problem()
# 	problem = Problem()
#
#     ### This should be done in one macro
# 	@variable(problem.model, grain[1:10] >= 0)
#     problem.vars[:grain] = (:control,size(grain))
#     ###
#
# 	@constraint(problem.model, sum(grain[i] for i in 1:10) <= 5)
#
#     ### This should be done in one macro
# 	@variable(problem.model, grain_prices[1:10])
#     problem.vars[:grain_prices] = (:random,size(grain_prices))
#     delete(problem.model,grain_prices)
#     ###
#
# 	# objective function is defined through the cutting plane oracle
# 	function my_cutting_plane_oracle()
#         function value(vars)::Float64
#
#         end
#         function gradient(vars)::Dict
#
#         end
#         return (value,gradient)
# 	end
#
# 	# distribution of random variables is defined through the measure oracle
# 	function my_measure_oracle()
#         function zero_order(vars_1,vars_1)::Float64
#         end
#         function first_order(vars_1,vars_2)::Dict
#         end
#         return (zero_order,first_order)
# 	end
#
#     problem.c_oracle = my_cutting_plane_oracle
#     problem.m_oracle = my_measure_oracle
#
# 	return problem
# end
#
# solve(BauckeAlgorithm(),build_problem(),Criteria(reltol))
#
# end
end
