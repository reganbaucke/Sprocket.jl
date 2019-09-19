####
# This file defines point and several related functions
####

# A point is a cartesian product of a point in control space and a point in the random space
struct Point
	control
	random
end

# returns the unit point for a given problem
function unit(problem)
	function new_ones(input)
		if input==()
			return 1.0
		else
			return ones(input)
		end
	end

	return map_over_keys(new_ones,problem[1]())
end

# returns the unit point for a given variable
function unit(problem,sym)
	base = origin(problem)
	point = problem[1]()

	if !isempty(point.control)
		for key in keys(point.control)
			if key == sym
				if !(typeof(point.control[key]) <: AbstractArray)
					base.control[key] = 1.0
					return base
				else
					base.control[key] = ones(point.control[key])
					return base
				end
			end
		end
	end

	if !isempty(point.random)
		for key in keys(point.random)
			if key == sym
				if !(typeof(point.random[key]) <: AbstractArray)
					base.random[key] = 1.0
					return base
				else
					base.random[key] = ones(point.random[key])
					return base
				end
			end
		end
	end
end

function unit(problem,sym,index)
	base = origin(problem)

	if !isempty(point.control)
		for key in keys(point.control)
			if key == sym
				if !(typeof(point.control[key]) <: AbstractArray)
					base.control[key] = 1.0
					return base
				else
					base.control[key][index] = 1.0
					return base
				end
			end
		end
	end

	if !isempty(point.random)
		for key in keys(point.random)
			if key == sym
				if !(typeof(point.random[key]) <: AbstractArray)
					base.random[key] = 1.0
					return base
				else
					base.random[key][index] = 1.0
					return base
				end
			end
		end
	end
end

# returns the origin point for a given problem
function origin(problem)
	function new_zeros(input)
		if input==()
			return 0.0
		else
			return zeros(input)
		end
	end

	return map_over_keys(new_zeros,problem[1]())
end


# define scalar multiplication
function Base.:*(scalar::Real,point::Point)
	function curried_multiply(arg)
		return scalar*arg
	end
	return map_over_point(curried_multiply, point)
end

## symmetry
function Base.:*(point::Point,scalar::Real)
	function curried_multiply(arg)
		return scalar*arg
	end
	return map_over_point(curried_multiply, point)
end

## define scalar divide
function Base.:/(point::Point,scalar::Real)
	return (1/scalar)*point
end

## define addition of points
function Base.:+(point_1::Point,point_2::Point)
	control = ()
	if !isempty(point_1.control)
		control = Dict()
		for key in keys(point_1.control)
			control[key] = point_1.control[key] + point_2.control[key]
		end
	end
	random = ()
	if !isempty(point_1.random)
		random = Dict()
		for key in keys(point_1.random)
			random[key] = point_1.random[key] + point_2.random[key]
		end
	end
	return Point(control,random)
end

## define subtractrion of points
function Base.:-(point_1::Point,point_2::Point)
	control = ()
	if !isempty(point_1.control)
		control = Dict()
		for key in keys(point_1.control)
			control[key] = point_1.control[key] - point_2.control[key]
		end
	end

	random = ()
	if !isempty(point_1.random)
		random = Dict()
		for key in keys(point_1.random)
			random[key] = point_1.random[key] - point_2.random[key]
		end
	end
	return Point(control,random)
end



#helper function, this will save a lot of copying and pasting
#return a new point with f function mapped over the elements of the current point
function map_over_point(f,point)
	newcontrol = ()
	#check control exists
	if !isempty(point.control)
		newcontrol = Dict()
		for key in keys(point.control)
			newcontrol[key] = map(f,point.control[key])
		end
	end

	newrandom = ()
	#check control exists
	if !isempty(point.random)
		newrandom = Dict()
		for key in keys(point.random)
			newrandom[key] = map(f,point.random[key])
		end
	end
	return Point(newcontrol,newrandom)
end

function map_over_keys(f,point)
	newcontrol = ()
	#check control exists
	if !isempty(point.control)
		newcontrol = Dict()
		newcontrol[key] = f(point.control[key])
	end

	newrandom = ()
	#check control exists
	if !isempty(point.random)
		newrandom = Dict()
		for key in keys(point.random)
			newrandom[key] = f(point.random[key])
		end
	end
	return Point(newcontrol,newrandom)
end
