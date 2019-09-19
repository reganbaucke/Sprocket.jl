struct Point
	control
	random
end

struct Cut
	point::Point
	value::Number
	grad::Point
end

struct Cut
	point::Point
	value::Number
	grad::Point
end

struct Problem
	type::Function
	m_oracle::Function
	c_oracle::Function
end

function solve(alg::Algorithm,prob::Problem,crit::Criteria)
	states{typeof(alg)}[]
	push!(states,initialise(alg,prob))
	while hasmet(alg::crit)
		push!(states,iterate(alg,state[end]))
	end
	return states
end

# best way is to input into solve all the functions that are called within solve.
# do not worry about dispatching on singleton types and calling functions outside
