####
# In this file we define the BauckeAlgorithm which will solve the optimisation problem exactly
####

using Combinatorics
# include("./Sprocket.jl")

module Baucke

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
		atom.corner_points = prob.domain

		println(atom.corner_points)

		atom.P = compute_probability(atom,prob.m_oracle)
		atom.A = compute_average_point(atom,prob.m_oracle)

		compute_weights!(atom)

		push!(atoms,atom)

		lower = Sprocket.Lowerbound(prob.vars,-99.0)
		upper = Sprocket.Upperbound(prob.vars,99.0,5.0)

		return (lower=lower,upper=upper,atoms=atoms,control=())
	end
	function iterate(state,prob)
		biggest_bound = largest_bound_gap(state.atoms,state.control,(state.lower,state.upper))
		center_point = compute_average_point(biggest_bound,prob.m_oracle)

		# new_state = deepcopy(state)
		atoms = copy(state.atoms)

		new_atoms = split_atom(biggest_bound,center_point)
		for atom in new_atoms
			atom.P = compute_probability(atom,prob.m_oracle)
			atom.A = compute_average_point(atom,prob.m_oracle)
			compute_weights!(atom)
		end

		union!(atoms,new_atoms)
		delete!(atoms,biggest_bound)

		cut = Sprocket.generate_cut(prob.c_oracle,center_point)

		Sprocket.update!(state.lower,cut)
		Sprocket.update!(state.upper,cut)

		# update_control_problem!(control_problem,cut)
		# get_new_control(control_problem)
		return (lower=state.lower,upper=state.upper,atoms=atoms,control=())
		# return new_state
	end
	function hasmet(crit::Sprocket.Criteria,state)
	end
	return (initialise,iterate,hasmet)
end

function compute_probability(atom::Baucke.Atom,oracle)
	return oracle(get_generating_pair(atom)...)[1]
end

function compute_average_point(atom::Baucke.Atom,oracle)
	return oracle(get_generating_pair(atom)...)[2]/oracle(get_generating_pair(atom)...)[1]
end

function compute_weights!(atom::Baucke.Atom)
	@assert is_bounded(atom)

	# A = zeros(dimension(vars)+1,length(atom.corner_points))
	A = Array{Float64}(undef,Sprocket.dimension(atom.A)+1,length(atom.corner_points))

	# save the indices in one place so the ordering over dicts with the keys is consistent
	indices = eachindex(atom.A)

	for (i,point) in enumerate(atom.corner_points)
		col = Float64[]
		for i in indices
			push!(col,point[i])
		end
		# push the last entry as the constraint that weights sum to one
		push!(col,1.0)
		A[:,i] = col
	end

	b = Float64[]
	for i in indices
		push!(b,atom.A[i])
	end
	# push the last entry as the constraint that weights sum to one
	push!(b,1.0)

	#compute the minimum norm whose solution corresponds to the most even weights
	weights = transpose(A)*(A*transpose(A)\b)


	atom.corner_weights = Dict{Any,Float64}()
	# assign weights
	for i in eachindex(atom.corner_points)
		atom.corner_weights[atom.corner_points[i]] = weights[i]
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
	upper_val = zeros(size(atoms))
	lower_val = zeros(size(atoms))
	# println("BOUND GAPS---------------")
	for (i,atom) in enumerate(atoms)
		upper_val[i] = upper_bound(atom,upper,control)
		lower_val[i] = lower_bound(atom,lower,control)
		gap[i] = upper_bound(atom,upper,control)-lower_bound(atom,lower,control)
		# println("--")
		# println(atom)
		# println(upper_val[i])
		# println(lower_val[i])

	end

	println("-----------")
	println(sum(gap))
	println("-")
	println(sum(upper_val))
	println(sum(lower_val))

	# println(lower.cuts)

	(value,index) = findmax(gap)

	# for (i,atom) in enumerate(atoms)
	# 	println("--")
	# 	if index == i
	# 		println("!!!!!")
	# 	end
	# 	println(atoms[i])
	# 	println(upper_val[i])
	# 	println(lower_val[i])
	# 	println("--")
	# end

	return atoms[index]
end

function lower_bound(atom,lower,control)
	return atom.P*(Sprocket.evaluate(lower,atom.A))
end

function is_bounded(atom::Baucke.Atom)
	for point in atom.corner_points
		if !all(map(x-> x!=Inf || x!=-Inf,point))
			return false
		end
	end
	return true
end


function upper_bound(atom,upper,control)
	# unbounded corner points?
	if !is_bounded(atom)
		return atom.P*(evaluate(upper,bounded_point(atom))+lipschitz*dist(bounded_point(atom),atom.A))
	end

	# bounded corner points
	sum = 0.0
	for point in atom.corner_points
		sum+= atom.corner_weights[point]*Sprocket.evaluate(upper,point)
	end

	return sum*atom.P
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

function get_generating_pair(a::Baucke.Atom)
	temp = sort(a.corner_points)
	return (temp[1],temp[end])
end


function rect_hull(x::Vector,y::Vector)
	hull=[]
	for i in powerset(1:length(x))
		z = copy(x)
		z[i] = y[i]
		push!(hull,z)
	end
	return hull
end

function rect_hull(x::Sprocket.Point,y::Sprocket.Point)
	# save in one place
	out = Sprocket.Point[]

	indices = eachindex(x)

	x_array = []
	y_array = []

	for i in indices
		push!(x_array,x[i])
		push!(y_array,y[i])
	end

	hull = rect_hull(x_array,y_array)

	for point in hull
		z = deepcopy(x)
		for i in 1:length(indices)
			z[indices[i]] = point[i]
		end
		push!(out,z)
	end
	return out
end


function Base.:<=(point_1::Sprocket.Point,point_2::Sprocket.Point)
	all(Sprocket.combine((x,y) -> x <= y, point_1,point_2))
end

function Base.:<(point_1::Sprocket.Point,point_2::Sprocket.Point)
	all(Sprocket.combine((x,y) -> x < y, point_1,point_2))
end

###
# for sorting
###
function Base.isless(point_1::Sprocket.Point,point_2::Sprocket.Point)
	point_1 <= point_2
end

function Base.all(point::Sprocket.Point)
	out = true
	for i in eachindex(point)
		out = out*point[i]
	end
	return out
end


function Base.string(atom::Baucke.Atom)
	out = ""
	for point in atom.corner_points
		out *= string(point) * "\n"
	end

	out*= "probability: $(atom.P) \n"
	out*= "barycenter: $(atom.A)"
end

function Base.show(io::IO,atom::Baucke.Atom)
	out = ""

	for point in atom.corner_points
		out *= string(point) * "\n"
	end

	out*= "probability: $(atom.P) \n"
	out*= "barycenter: $(atom.A)"
	println(io,out)
end

# TODO
# add way of computing upper bound for unbounded atom
