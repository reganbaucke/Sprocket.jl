include("../src/Sprocket.jl")
include("../src/baucke.jl")
# using Sprocket
using JuMP

##
# One dimensional problem, no control
##

function test_one()

function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	xi = Sprocket.Variable(name=:xi, size=(),type=Sprocket.Random())
	Sprocket.add_variable(problem,xi)

	#set the domain of the random variable
	Sprocket.set_domain(problem,rect_hull,(Sprocket.Point(Dict(xi =>0.0)), Sprocket.Point(Dict(xi => 1.0))))

	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function oracle(vars)
			value = vars[:xi]^4
			grad = deepcopy(vars)
			grad[:xi] = 4*vars[:xi]^3
			return (value,grad)
		end
		return oracle
	end

	# distribution of random variables is defined through the measure oracle
	function my_measure_oracle()
		function oracle(vars_1,vars_2)
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
			fresh = deepcopy(vars_1)
			fresh[:xi] = (right^2 - left^2)/2
			return (right-left,fresh)
		end
		return oracle
	end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

	return problem
end


my_prob = build_problem()
Sprocket.Criteria(reltol=0.1)
states = Sprocket.solve(BauckeAlgorithm(),my_prob,Sprocket.Criteria(reltol=0.1))

return nothing
### ANSWER = 0.2
end



##
# Test 2, two dimensional problem no control
##

function test_two()

function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	xi = Sprocket.Variable(name=:xi, size=(),type = Sprocket.Random())
	Sprocket.add_variable(problem,xi)

	# add another
	lam = Sprocket.Variable(name=:lam, size=(),type = Sprocket.Random())
	Sprocket.add_variable(problem,lam)

	#set the domain of the random variable
	Sprocket.set_domain(problem,rect_hull,(Sprocket.Point(Dict(xi =>0.0,lam => 0.0 )), Sprocket.Point(Dict(xi => 1.0, lam => 1.0))))

	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function oracle(vars)
			value = vars[:xi]^2 + 3*vars[:lam]^2 + 0.1*vars[:lam]*vars[:xi]

			grad = deepcopy(vars)
			grad[:xi] = 2*vars[:xi] + 0.1*vars[:lam]
			grad[:lam] = 0.1*vars[:xi] + 6*vars[:lam]

			return (value,grad)
		end
		return oracle
	end

	# distribution of random variables is defined through the measure oracle
	function my_measure_oracle()
		function oracle(vars_1,vars_2)
			left_xi = vars_1[:xi]
			right_xi = vars_2[:xi]
			if vars_1[:xi] <= 0
				left_xi = 0
			end
			if vars_1[:xi] >= 1
				left_xi = 1
			end
			if vars_2[:xi] <= 0
				right_xi = 0
			end
			if vars_2[:xi] >= 1
				right_xi = 1
			end
			factor_xi = right_xi - left_xi

			left_lam = vars_1[:lam]
			right_lam = vars_2[:lam]
			if vars_1[:lam] <= 0
				left_lam = 0
			end
			if vars_1[:lam] >= 1
				left_lam = 1
			end
			if vars_2[:lam] <= 0
				right_lam = 0
			end
			if vars_2[:lam] >= 1
				right_lam = 1
			end
			factor_lam = right_lam - left_lam

			fresh = deepcopy(vars_1)
			fresh[:xi] = (right_xi^2 - left_xi^2)*factor_lam/2
			fresh[:lam] = (right_lam^2 - left_lam^2)*factor_xi/2
			return (factor_lam*factor_xi,fresh)
		end
		return oracle
	end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

	return problem
end

my_prob = build_problem()
# Sprocket.Criteria(reltol=0.1)
states = Sprocket.solve(BauckeAlgorithm(),my_prob,Sprocket.Criteria(reltol=0.1))

# ANSWER = 1.35833
end

function test_three()
function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	xi = Sprocket.Variable(name=:xi, size=(),type = Sprocket.Random())
	Sprocket.add_variable(problem,xi)

	# add another variable to the problem, this one is an 1 dimensional array of size 2
	lam = Sprocket.Variable(name=:lam, size=(2), type = Sprocket.Random())
	Sprocket.add_variable(problem,lam)

	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function oracle(vars)
			value = vars[:xi]^2 + 3*vars[:lam]^2 + 0.1*vars[:lam]*vars[:xi]

			grad = deepcopy(vars)
			grad[:xi] = 2*vars[:xi] + 0.1*vars[:lam]
			grad[:lam] = 0.1*vars[:xi] + 6*vars[:lam]

			return (value,grad)
		end
		return oracle
	end

	# distribution of random variables is defined through the measure oracle
	function my_measure_oracle()
		function oracle(vars_1,vars_2)
			left_xi = vars_1[:xi]
			right_xi = vars_2[:xi]
			if vars_1[:xi] <= 0
				left_xi = 0
			end
			if vars_1[:xi] >= 1
				left_xi = 1
			end
			if vars_2[:xi] <= 0
				right_xi = 0
			end
			if vars_2[:xi] >= 1
				right_xi = 1
			end
			factor_xi = right_xi - left_xi

			left_lam = vars_1[:lam]
			right_lam = vars_2[:lam]
			if vars_1[:lam] <= 0
				left_lam = 0
			end
			if vars_1[:lam] >= 1
				left_lam = 1
			end
			if vars_2[:lam] <= 0
				right_lam = 0
			end
			if vars_2[:lam] >= 1
				right_lam = 1
			end
			factor_lam = right_lam - left_lam

			fresh = deepcopy(vars_1)
			fresh[:xi] = (right_xi^2 - left_xi^2)*factor_lam/2
			fresh[:lam] = (right_lam^2 - left_lam^2)*factor_xi/2
			return (factor_lam*factor_xi,fresh)
		end
		return oracle
	end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

	return problem
end

my_prob = build_problem()
# Sprocket.Criteria(reltol=0.1)
states = Sprocket.solve(BauckeAlgorithm(),my_prob,Sprocket.Criteria(reltol=0.1))

# ANSWER = 1.35833
end
