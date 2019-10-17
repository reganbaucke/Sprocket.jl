include("../src/Sprocket.jl")
include("../src/Baucke.jl")
using .Sprocket
using .Baucke
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
	Sprocket.set_domain(problem,Baucke.rect_hull,(Sprocket.Point(Dict(xi =>0.0)), Sprocket.Point(Dict(xi => 1.0))))

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
states = Sprocket.solve(Baucke.BauckeAlgorithm(),my_prob,Sprocket.Criteria(iterations=10))

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
	Sprocket.set_domain(problem,Baucke.rect_hull,(Sprocket.Point(Dict(xi =>0.0,lam => 0.0 )), Sprocket.Point(Dict(xi => 1.0, lam => 1.0))))

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
states = Sprocket.solve(Baucke.BauckeAlgorithm(),my_prob,Sprocket.Criteria(iterations = 10))
# ANSWER = 1.35833
end

##
# Test 3, two dimensional problem, one variable with two entries
##

function test_three()
function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	xi = Sprocket.Variable(name=:xi, size=(2,) ,type = Sprocket.Random())
	Sprocket.add_variable(problem,xi)

	#set the domain of the random variable
	Sprocket.set_domain(problem,Baucke.rect_hull,(Sprocket.Point(Dict(xi =>[0.0,0.0])), Sprocket.Point(Dict(xi => [1.0,1.0]))))

	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function oracle(vars)
			value = vars[:xi][1]^2 + 3*vars[:xi][2]^2 + 0.1*vars[:xi][1]*vars[:xi][2]

			grad = deepcopy(vars)
			grad[:xi][1] =   2*vars[:xi][1] + 0.1*vars[:xi][2]
			grad[:xi][2] = 0.1*vars[:xi][1] +   6*vars[:xi][2]

			return (value,grad)
		end
		return oracle
	end

	# distribution of random variables is defined through the measure oracle
	function my_measure_oracle()
		function oracle(vars_1,vars_2)
			left_xi = vars_1[:xi][1]
			right_xi = vars_2[:xi][1]
			if vars_1[:xi][1] <= 0
				left_xi = 0
			end
			if vars_1[:xi][1] >= 1
				left_xi = 1
			end
			if vars_2[:xi][1] <= 0
				right_xi = 0
			end
			if vars_2[:xi][1] >= 1
				right_xi = 1
			end
			factor_xi = right_xi - left_xi

			left_lam = vars_1[:xi][2]
			right_lam = vars_2[:xi][2]
			if vars_1[:xi][2] <= 0
				left_lam = 0
			end
			if vars_1[:xi][2] >= 1
				left_lam = 1
			end
			if vars_2[:xi][2] <= 0
				right_lam = 0
			end
			if vars_2[:xi][2] >= 1
				right_lam = 1
			end
			factor_lam = right_lam - left_lam

			fresh = deepcopy(vars_1)
			fresh[:xi][1] = (right_xi^2 - left_xi^2)*factor_lam/2
			fresh[:xi][2] = (right_lam^2 - left_lam^2)*factor_xi/2
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
states = Sprocket.solve(Baucke.BauckeAlgorithm(),my_prob,Sprocket.Criteria(iterations = 60))

return nothing

# ANSWER = 1.35833
end

##
# One random variable, one control
##
function test_four()

function build_problem()
	#create an empty problem
	problem = Sprocket.Problem()

	# add one real valued random variable to the problem
	xi = Sprocket.Variable(name=:xi, size=(),type=Sprocket.Random())
	Sprocket.add_variable(problem,xi)

	#set the domain of the random variable
	Sprocket.set_domain(problem,Baucke.rect_hull,(Sprocket.Point(Dict(xi =>0.0)), Sprocket.Point(Dict(xi => 1.0))))

	u = Sprocket.Variable(name=:u, size=(),type=Sprocket.Control())
	add_variable(problem,u)
	@constraint(problem.model, -1.0 <= u <= 1.0)


	# objective function is defined through the cutting plane oracle
	function my_cutting_plane_oracle()
		function oracle(vars)
			value = vars[:xi]^2 + 0.1*vars[:u]*vars[:xi] + vars[:u]^2
			grad = deepcopy(vars)
			grad[:xi] = 2.0*vars[:xi] + 0.1*vars[:u]
			grad[:u]  = 0.1*vars[:xi] + 2.0*vars[:u]
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
states = Sprocket.solve(Baucke.BauckeAlgorithm(),my_prob,Sprocket.Criteria(iterations = 10))

return nothing
### ANSWER = 0.2
end


function test_all()
	test_one()
	test_two()
	test_three()
	# test_four()
end
