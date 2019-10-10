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
	var_dict::Dict{Variable, Any}
end

function Point(dict::Dict{T,Any}) where T <: Variable
	my_dict = Dict{Variable,Any}()
	for key in keys(dict)
		my_dict = copy(dict[key])
	end
	return Point(my_dict)
end

function Base.eachindex(point::Point)
	vars = collect(keys(point.var_dict))
	names = map(x->x.name,vars)
	out = []
	for i in 1:length(vars)
		if size(point.var_dict[vars[i]]) == ()
			push!(out,(names[i],))
		else
			for j = 1:length(point.var_dict[vars[i]])
				push!(out,(names[i],j))
			end
		end
	end
	return out
end

function Base.eachindex(vars::Set{Variable})
	vars = collect(vars)
	names = map(x->x.name,vars)
	out = []
	for i in 1:length(vars)
		if vars[i].size == ()
			push!(out,(names[i],))
		else
			for j = 1:prod(collect(vars[i].size))
				push!(out,(names[i],j))
			end
		end
	end
	return out
end

function Base.getindex(my_point::Point,i)
	if typeof(i) == Symbol
		return my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i]]
		# filter(x-> x.name==i[1],collect(keys(my_point.var_dict)))[1]
	end
	if length(i) == 1
		return my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i[1]]]
	end
	if length(i) == 2
		return my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i[1]]][i[2]]
	end
	# return my_point.var_dict[my_var]
end

function Base.setindex!(my_point::Point,data,i)
	if typeof(i) == Symbol
		my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i]] = data
		return nothing
	end
	if length(i) == 1
		my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i[1]]] = data
		return nothing
	end
	if length(i) == 2
		my_point.var_dict[Set(collect(keys(my_point.var_dict)))[i[1]]][i[2]] = data
		return nothing
	end
	# return my_point.var_dict[my_var]
end

function Base.getindex(vars::Set{Variable},i)
	# @assert length(i) == 1
	return filter(x-> x.name==i,collect(vars))[1]
	# return my_point.var_dict[my_var]
end

function get_vars(point::Point)
	return Set(collect(keys(point.var_dict)))
end

function Base.iterate(point::Point)
	if isempty(point.var_dict)
		return nothing
	end
	vars = collect(keys(point.var_dict))
	if typeof(point.var_dict[vars[1]]) <: AbstractArray
		return (point.var_dict[vars[1]][1],(vars=vars,prev_lin_ind = 1))
	else
		return (point.var_dict[vars[1]],(vars=vars[2:end],prev_lin_ind=0))
	end
end

function Base.iterate(point::Point,state)
	if isempty(state.vars)
		return nothing
	end
	if typeof(point.var_dict[state.vars[1]]) <: AbstractArray
		if (state.prev_lin_ind + 1) == length(point.var_dict[state.vars[1]])
			return (point.var_dict[state.vars[1]][end],(vars = state.vars[2:end],prev_lin_ind = 0))
		else
			return (point.var_dict[state.vars[1]][state.prev_lin_ind + 1],(vars=state.vars,prev_lin_ind = state.prev_lin_ind + 1))
		end
	else
		return (point.var_dict[state.vars[1]],(vars = state.vars[2:end],prev_lin_ind=0))
	end
end

function Base.map(f,point::Point)
	out = deepcopy(point)
	for var in keys(point.var_dict)
		out.var_dict[var] = f.(point.var_dict[var])
		if typeof(var) <: AbstractArray
			for (val,i) in enumerate(point.var_dict[var])
				out.var_dict[var][i] = f(val)
			end
		else
			out.var_dict[var] = f(point.var_dict[var])
		end
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

function Base.:*(var_1::Variable,var_2::Variable)
	return Set([var_1,var_2])
end

function Base.:*(var_1::Set{Variable},var_2::Variable)
	return union(var_1,var_2)
end

function Base.:*(var_1::Variable,var_2::Set{Variable})
	# check that none of the keys are the same
	return union(var_2,var_1)
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

hello = Variable(:xi,(2,),Control())
blue = Variable(:melt,(),Random())

ultra = Point(Dict(hello => [10.0,80.0]))
thick = Point(Dict(blue => 9.0))

my_point = Point(Dict(hello => 9.0))
new_point = map(x->x+1,my_point)
another_point = combine((x,y) -> x + y, my_point, new_point)

my_point = Point(Dict(hello => 9.0))
