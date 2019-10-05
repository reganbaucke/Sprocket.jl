using Test
using JuMP
using GLPK
# using Sprocket
include("../src/Sprocket.jl")


function build_type(point)
	function type()
		return point
	end
end

function build_c_oracle()
	function zero_order(vars)
		return vars[:xi][2]^2
	end
	function first_order(vars)
		grad = Dict()
		grad[:xi] = (vars[:xi][1], 2*vars[:xi][2])
		return grad
	end
	return (zero_order,first_order)
end

function build_m_oracle()
	function zero_order(vars_1,vars_2)
		left = -Inf
		right = Inf
		if random_1[:xi] <= 0
			left = 0
		end
		if random_1[:xi] >= 1
			left = 1
		end
		if random_2[:xi] <= 0
			right = 0
		end
		if random_2[:xi] >= 1
			right = 1
		end
		return right-left
	end
	function first_order(vars_1,vars_2)
		left = -Inf
		right = Inf
		if random_1[:xi][2] <= 0
			left = 0
		end
		if random_1[:xi][2] >= 1
			left = 1
		end
		if random_2[:xi][2] <= 0
			right = 0
		end
		if random_2[:xi][2] >= 1
			right = 1
		end
		fresh = copy(random_1)
		fresh[:xi] = (right^2 - left^2)/2
		return fresh
	end
	return (zero_order,first_order)
end

# make a problem
my_prob = build_problem(build_point())


vars = Dict()

vars[:xi] = (:random,())

# make a lowerbound
my_lower = Lowerbound(vars,-999.9)
my_upper = Upperbound(vars,999.9,10)

origin = Dict()
origin[:xi] = (:random, 0.0)

unit = Dict()
unit[:xi] = (:random, 1.0)

evaluate(my_lower,origin)
evaluate(my_upper,origin)

my_cut = generate_cut(build_c_oracle(),unit)

update!(my_lower,my_cut)
update!(my_upper,my_cut)
