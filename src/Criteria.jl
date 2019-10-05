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
