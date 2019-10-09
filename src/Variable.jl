abstract type VariableType end
struct Control <: VariableType end
struct Random <: VariableType end

struct Variable{T <: VariableType}
	name::Symbol
	size::Tuple
	type::T
end

# Wrapper for dictionary
struct Point
	var_dict::Dict{Variable,Union{Float64, Array{Float64}}}
end

function Base.getindex(my_point::Point,my_var::Variable)
	return my_point.var_dict[my_var]
end

function Base.map(f,point::Point)
	out = deepcopy(point)
	for var in keys(point.var_dict)
		out.var_dict[var] = f.(point.var_dict[var])
	end
	return out
end

function combine(f,point_1::Point,point_2::Point)
	out = deepcopy(point_1)
	for var in keys(point_1.var_dict)
		out.var_dict[var] = f.(point_1.var_dict[var],point_2.var_dict[var])
	end
	return out
end

function fold(f,initial_value,point::Point)
	my_initial = copy(initial_value)
	for var in keys(point.var_dict)
		if typeof(var) <: AbstractArray
			for val in point.var_dict[var]
				my_initial = f(my_initial,val)
			end
		else
			my_initial = f(my_initial,point.var_dict[var])
		end
	end
	return my_initial
end


##
# Basic operations on points
##
function Base.:*(point_1::Point,point_2::Point)
	# check that none of the keys are the same
	out = deepcopy(point_1)
	for var in keys(point_2.var_dict)
		out.var_dict[var] = copy(point_2.var_dict[var])
	end
	return out
end

function Base.:+(point_1::Point,point_2::Point)
	combine(+, point_1,point_2)
end

function Base.:-(point_1::Point,point_2::Point)
	combine(-, point_1,point_2)
end

function Base.:-(point::Point)
	map(-, point)
end

function l_one_norm(my_point::Point)
	fold((x,y)->x + abs(y),0.0,my_point)
end

function Base.:*(my_point::Point,scalar::Number)
	map(x->x*scalar,my_point)
end

Base.:*(scalar::Number,my_point::Point) = Base.:*(my_point::Point,scalar::Number)

function Base.:/(my_point::Point,scalar::Number)
	map(x->x*(1/scalar),my_point)
end


## basic tests
# hello = d.ControlVariable(:xi,())
#
# my_point = d.Point(Dict(hello => 9.0))
# new_point = map(x->x+1,my_point)
# another_point = d.combine((x,y) -> x + y, my_point, new_point)
#
# my_point = d.Point(Dict(hello => 9.0))


## basic tests

# hello = d.Variable(:xi,(),d.Control())
# blue = d.Variable(:melt,(),d.Random())

# ultra = d.Point(Dict(hello => 10.0))
# thick = d.Point(Dict(blue => 9.0))

# my_point = d.Point(Dict(hello => 9.0))
# new_point = map(x->x+1,my_point)
# another_point = d.combine((x,y) -> x + y, my_point, new_point)

# my_point = d.Point(Dict(hello => 9.0))
