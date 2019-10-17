####
# This file defines the criteria type and related functions
####

struct Criteria
	# tolerances
	abstol::Real
	reltol::Real

	# time
	time::Real

	#iterations
	iterations::Int64

	# number of calls to the oracles
	m_calls::Int64
	c_calls::Int64
end

### Constructor with defaults
function Criteria(;abstol::Float64=0.0,
	 			   reltol::Float64=0.0,
				   time::Float64=Inf,
				   iterations::Int64=9223372036854775807,
				   m_calls::Int64=9223372036854775807,
				   c_calls::Int64=9223372036854775807)
	 return Criteria(abstol,reltol,time,iterations,m_calls,c_calls)
end

function initial_criteria()
	 return Criteria(abstol = Inf,reltol = Inf,time = 0.0,iterations = 0,m_calls = 0,c_calls =0)
end

function met(crit::Criteria, target::Criteria)
	if crit.abstol <= target.abstol
		return true
	end

	if crit.reltol <= target.reltol
		return true
	end

	if crit.time >= target.time
		return true
	end

	if crit.iterations >= target.iterations
		return true
	end

	if crit.m_calls >= target.m_calls
		return true
	end

	if crit.c_calls >= target.c_calls
		return true
	end
	return false
end

function id(x)
	x
end

function update_with(old_crit;abstol=id, reltol=id, time=id, iterations=id, m_calls=id, c_calls=id)
	return Criteria(
		abstol(old_crit.abstol),
		reltol(old_crit.reltol),
		time(old_crit.time),
		iterations(old_crit.iterations),
		m_calls(old_crit.m_calls),
		c_calls(old_crit.c_calls))
end
