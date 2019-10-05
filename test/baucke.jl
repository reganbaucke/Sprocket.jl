include("../src/Sprocket.jl")
include("../src/baucke.jl")
# using Sprocket
using JuMP

function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	@variable(problem.model, xi)
    problem.vars[:xi] = (:random,size(xi))
    delete(problem.model,xi)
    ###

	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function zero_order(vars)
			return vars[:xi]^4
			# return 0
		end
		function first_order(vars)
			grad = Dict()
			grad[:xi] = 4*vars[:xi]^3
			# grad[:xi] = 0
			return grad
		end
		return (zero_order,first_order)
	end

	# distribution of random variables is defined through the measure oracle
	function my_measure_oracle()
		function zero_order(vars_1,vars_2)
			left = vars_1[:xi]
			right = vars_2[:xi]
			if vars_1[:xi] <= 0
				left = 0
			end
			if vars_1[:xi] >= 1
				left = 1
			end
			if vars_2[:xi] <= 0
				right = 0
			end
			if vars_2[:xi] >= 1
				right = 1
			end
			return right-left
		end
		function first_order(vars_1,vars_2)
			left = vars_1[:xi]
			right = vars_2[:xi]
			if vars_1[:xi] <= 0
				left = 0
			end
			if vars_1[:xi] >= 1
				left = 1
			end
			if vars_2[:xi] <= 0
				right = 0
			end
			if vars_2[:xi] >= 1
				right = 1
			end
			fresh = copy(vars_1)
			fresh[:xi] = (right^2 - left^2)/2
			return fresh
		end
		return (zero_order,first_order)
	end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

	return problem
end

my_prob = build_problem()

# Sprocket.Criteria(reltol=0.1)

states = Sprocket.solve(BauckeAlgorithm(),my_prob,Sprocket.Criteria(reltol=0.1))
