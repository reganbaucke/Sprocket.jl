using Test
using JuMP
using GLPK
using Sprocket

Base.copy(()) = ()

function build_point()
	control = ()

	random = Dict()

	# singleton
	random[:xi] = ()

	# one dimensional array of length 5
	return Point(control,random)
end

function build_problem(point::Point)
	cutting_oracle = build_c_oracle()
	m_oracle = build_m_oracle()
	type = build_type(point)
	return (type,cutting_oracle,m_oracle)
end

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
	function zero_order(random_1,random_2)
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
	function first_order(random_1,random_2)
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
