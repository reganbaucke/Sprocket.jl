####
# In this file we define the BauckeAlgorithm which will solve the optimisation problem exactly
####

using Combinatorics

module Baucke
# using Sprocket
include("./Sprocket.jl")

#structs for use in BauckeAlgorithm
mutable struct Atom
	corner_points
	corner_weights
	P
	A
	function Atom()
		return new()
	end
end

end


function BauckeAlgorithm()
	function initialise(prob)

		atoms = Set()

		atom = Baucke.Atom()
		atom.corner_points = []
		point_1 = Dict()
		point_1[:xi] = 0.0
		point_2 = Dict()
		point_2[:xi] = 1.0
		push!(atom.corner_points,point_1)
		push!(atom.corner_points,point_2)

		atom.corner_weights = Dict()
		atom.corner_weights[point_1] = 0.5
		atom.corner_weights[point_2] = 0.5

		atom.A = Dict(:xi => 0.5)
		atom.P = 1.0

		push!(atoms,atom)

		lower = Sprocket.Lowerbound(prob.vars,-99)
		upper = Sprocket.Upperbound(prob.vars,99,10)
		return (lower=lower,upper=upper,atoms=atoms,control=())
	end
	function iterate(state,prob)
		println(state.lower)

		biggest_bound = largest_bound_gap(state.atoms,state.control,(state.lower,state.upper))
		center_point = compute_average_point(biggest_bound,prob.m_oracle)
		atoms = state.atoms

		union!(atoms,split_atom(biggest_bound,center_point))
		delete!(atoms,biggest_bound)

		cut = Sprocket.generate_cut(prob.c_oracle,center_point)
		Sprocket.update!(state.lower,cut)
		Sprocket.update!(state.upper,cut)

		# update_control_problem!(control_problem,cut)
		# get_new_control(control_problem)
		return (lower=state.lower,upper=state.upper,atoms=state.atoms,control=())
	end
	function hasmet(crit::Sprocket.Criteria,state)
	end
	return (initialise,iterate,hasmet)
end

function compute_probability(atom::Baucke.Atom,oracle)
	return m_oracle[1](get_generating_pair(atom)...)
end

function compute_average_point(atom::Baucke.Atom,oracle)
	return oracle[2](get_generating_pair(atom)...)/oracle[1](get_generating_pair(atom)...)
end

function compute_weights(atom::Baucke.Atom)
	@assert is_bounded(atom)

	# save the keys in one place so the ordering over dicts with the keys is consistent
	var_keys = keys(atom.A)
	# A = zeros(dimension(vars)+1,length(atom.corner_points))
	A = Array{Float64}(undef,dimension(vars)+1,length(atom.corner_points))

	for (point,i) in enumerate(atom.corner_points)
		col = Float64[]
		for key in var_keys
			if point[key] <: Number
				push!(col,point[key])
			else
				for j in eachindex(point[key])
					push!(col,point[key][j])
				end
			end
		end
		# push the last entry as the constraint that weights sum to one
		push!(col,1.0)
		A[:,i] = col
	end

	b = Float64[]
	for key in var_keys
		if point[key] <: Number
			push!(b,atom.A[key])
		else
			for j in eachindex(atom.A[key])
				push!(col,atom.A[key][j])
			end
		end
	end

	# push the last entry as the constraint that weights sum to one
	push!(col,1.0)

	#compute the minimum norm whose solution corresponds to the most even weights
	weights = transpose(A)*(A*transpose(A)\b)

	# assign weights
	for i in eachindex(atom.corner_points)
		atom.compute_weights[atom.corner_points] = weights[i]
	end

	return nothing
end

function split_atom(atom::Baucke.Atom,center)
	atoms = Set()
	for corner in atom.corner_points
		new_atom = Baucke.Atom()
		new_atom.corner_points = rect_hull(corner,center)
		push!(atoms,new_atom)
	end
	return atoms
end


function largest_bound_gap(atoms,control,(lower,upper))
	atoms = collect(atoms)
	gap = zeros(size(atoms))
	for (i,atom) in enumerate(atoms)
		gap[i] = upper_bound(atom,upper,control)-lower_bound(atom,lower,control)
	end
	(value,index) = findmax(gap)
	return atoms[index]
end

function lower_bound(atom,lower,control)
	return atom.P*(Sprocket.evaluate(lower,atom.A))
end

function is_bounded(atom::Baucke.Atom)
	for point in atom.corner_points
		if !all(map_over_vars(x-> x!=Inf || x!=-Inf,point))
			return false
		end
	end
	return true
end


function upper_bound(atom,upper,control)
	# unbounded corner points?
	if !is_bounded(atom)
		return atom.P*(evaluate(upper,control*bounded_point(atom))+lipschitz*dist(bounded_point(atom),atom.A))
	end

	# bounded corner points
	sum = 0.0
	for point in atom.corner_points
		sum+= atom.corner_weights[point]*Sprocket.evaluate(upper,point)
	end

	return sum
end

function get_probability(m_oracle,a::Baucke.Atom)
	(zero,first) = m_oracle
	return zero(get_generating_pair(a...))/zero(get_generating_pair(a...))
end

function get_probability(m_oracle,a::Baucke.Atom)
	(zero,first) = m_oracle
	return zero(get_generating_pair(a...))/zero(get_generating_pair(a...))
end

function get_average_point(m_oracle,a::Baucke.Atom)
	(zero,first) = m_oracle
	return first(get_generating_pair(a...))
end

function compute_lower_bound(lower,atom,vars)

end

function compute_upper_bound(lower,atom,vars)

end

function compute_bound_gap(lower,atom,vars)

end

function get_generating_pair(a::Baucke.Atom)
	temp = sort(a.corner_points)
	return (temp[1],temp[end])
end


function map_over_vars(f,vars)
	output = Dict()

	for key in keys(vars)
		output[key] = map(f,vars[key])
	end

	return output
end


###
# This should be a symmetric function
###
function applic_over_vars(f,var_1,var_2)
	output = Dict()

	for key in keys(var_1)
		output[key] = f.(var_1[key],var_2[key])
	end

	return output
end

function fold_over_vars(f,init,var)
	output = copy(init)

	for key in keys(var)
		if !(typeof(var[key]) <: AbstractArray)
			output = f(output,var[key])
		else
			for element in var[key]
				output = f(output,element)
			end
		end
	end

	return output
end

function map_over_keys(f,vars)
	output = Dict()

	for key in keys(vars)
		output[key] = f(vars[key])
	end

	return output
end

# define scalar multiplication
function Base.:*(scalar::Real,vars::Dict)
	function curried_multiply(arg)
		return scalar*arg
	end
	return map_over_vars(curried_multiply, vars)
end



## symmetry
function Base.:*(vars::Dict,scalar::Real)
	Base.:*(scalar,vars)
end

## define scalar divide
function Base.:/(vars::Dict,scalar::Real)
	return (1/scalar)*vars
end

## define addition of points
function Base.:+(vars_1::Dict,vars_2::Dict)
	applic_over_vars(+,vars_1,vars_2)
end

## define subtractrion of points
function Base.:-(vars_1,vars_2)
	applic_over_vars(-,vars_1::Dict,vars_2::Dict)
end


function rect_hull(x,y)
	hull=[]
	for i in powerset(1:length(x))
		z = copy(x)
		z[i] = y[i]
		push!(hull,z)
	end
	return hull
end

function flatten(dict,keys)
	flat = Array{Tuple,1}()
	for key in keys
		if !(typeof(dict[key]) <: AbstractArray)
			push!(flat,(key,dict[key]))
		else
			for (i,v) in enumerate(dict[key])
				push!(flat,(key,i,size(dict[key]),dict[key][i]))
			end
		end
	end
	return flat
end

function dict_zip(input)
	my_dict = Dict()
	for each in input
		if length(each) == 2
			my_dict[each[1]] = each[2]
		else
			if !in(each[1],keys(my_dict))
				my_dict[each[1]] = NaN*ones(each[3])
				my_dict[each[1]][each[2]] = each[4]
			else
				my_dict[each[1]][each[2]] = each[4]
			end
		end
	end
	return my_dict
end

function rect_hull(x::Dict,y::Dict)
	the_keys = keys(x)
	flat_points = rect_hull(flatten(x,the_keys),flatten(y,the_keys))
	out = Set()
	for point in flat_points
		push!(out,dict_zip(point))
	end
	return collect(out)
end

function one_norm(x::Dict)
	return fold_over_vars((cum,element) -> cum + abs(element),0.0,x)
end

function Base.:<=(var_1::Dict,var_2::Dict)
	all(map_over_keys(all,applic_over_vars((x,y) -> all(x .<= y ),var_1,var_2)))
end

###
# for sorting
###
function Base.isless(var_1::Dict,var_2::Dict)
	var_1 <= var_2
end

function Base.all(vars::Dict)
	out = true
	for key in keys(vars)
		out = out*vars[key]
	end
	return out
end

function dimension(var::Dict)
	fold_over_vars((x,y)-> x+1,0,var)
end


function Base.:<(var_1::Dict,var_2::Dict)
	all(map_over_keys(all,applic_over_vars((x,y) -> all(x .< y ),var_1,var_2)))
end

# TODO
# add the corner weights function
