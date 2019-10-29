####
# In this file we define the BauckeAlgorithm which will solve the optimisation problem exactly
####
module Baucke

using Combinatorics
using ..Sprocket
using JuMP

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

include("./ControlProblem.jl")

struct State
	lower::Sprocket.Lowerbound
	upper::Sprocket.Upperbound
	atoms::Set{Baucke.Atom}
	control
	criteria::Sprocket.Criteria
end

function BauckeAlgorithm()
	function initialise(prob)

		atoms = Set()

		atom = Baucke.Atom()
		atom.corner_points = prob.domain

		atom.P = compute_probability(atom,prob.m_oracle)
		atom.A = compute_average_point(atom,prob.m_oracle)

    if atom |> is_bounded
  		compute_weights!(atom)
    end

		push!(atoms,atom)

		lower = Sprocket.Lowerbound(prob.vars,-99.0)
		upper = Sprocket.Upperbound(prob.vars,99.0,5.0)

		###
		# Set up the control problem
		###
		control_problem = build_control_problem(prob)
		if control_problem != nothing
			add_atom(control_problem.value , atom)
		end

		###
		# Set up the initial status of the criteria
		###
		new_criteria = Sprocket.initial_criteria()

		return Baucke.State(lower,upper,atoms,control_problem,new_criteria)
	end


	function iterate(state,prob)

		control = map(get_new_control, state.control)

		# compute the bound gap at the new control
		# biggest_bound = largest_bound_gap(state.atoms,(state.lower,state.upper))
		(biggest_bound, upper_bound, lower_bound )= largest_bound_gap(state.atoms,control,(state.lower,state.upper))
		center_point = compute_average_point(biggest_bound,prob.m_oracle)

		# new_state = deepcopy(state)
		atoms = copy(state.atoms)

		new_atoms = split_atom(biggest_bound,center_point)
		for atom in new_atoms
			atom.P = compute_probability(atom,prob.m_oracle)
			atom.A = compute_average_point(atom,prob.m_oracle)
      if atom |> is_bounded
  			compute_weights!(atom)
      end
		end

		# update collection of atoms
		union!(atoms,new_atoms)
		delete!(atoms,biggest_bound)

		# get cut at this point
		cut = Sprocket.generate_cut(prob.c_oracle,center_point*control)

		# update control problem
    if state.control != nothing
  		delete_atom(state.control.value,biggest_bound)
  		for atom in new_atoms
  			add_atom(state.control.value,atom)
  			gen_constraints_for_new_atom(state.control.value,atom)
  		end
  		update_all_atoms_with_cut(state.control.value,cut)
  		add_cut!(state.control.value,cut)
    end

		Sprocket.update!(state.lower,cut)
		Sprocket.update!(state.upper,cut)

		###
		# Update Criteria
		###

		# new_criteria = Sprocket.update_with(state.criteria,
    #   iterations = x -> x + 1)
    #   iterations = x->x+1)

    new_criteria = Sprocket.update_with(state.criteria,
      iterations = state.criteria.iterations + 1,
      upper = upper_bound,
      lower = lower_bound,
      abstol = abstol(upper_bound,lower_bound),
      reltol = reltol_upper(upper_bound,lower_bound)
    )

		# return new_state
		return Baucke.State(state.lower,state.upper,atoms,state.control,new_criteria)
	end

  function hasmet(crit::Sprocket.Criteria,state)
    Sprocket.met(state.criteria,crit)
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
	for (i,atom) in enumerate(atoms)
		upper_val[i] = upper_bound(atom,upper,control)
		lower_val[i] = lower_bound(atom,lower,control)
		gap[i] = upper_bound(atom,upper,control)-lower_bound(atom,lower,control)
	end

	(value,index) = findmax(gap)

	return (atoms[index], sum(upper_val), sum(lower_val))
end

function lower_bound(atom,lower,control)
		atom.P*(Sprocket.evaluate(lower,atom.A*control))
end

function is_bounded(atom::Baucke.Atom)
	for point in atom.corner_points
		if any(map(x -> x == Inf || x == -Inf, point))
			return false
		end
	end
	return true
end

function is_fully_unbounded(atom::Baucke.Atom)
	for point in get_generating_pair(atom)
		if !all(map(x -> x == Inf || x == -Inf, point))
			return false
		end
	end
	return true
end


function upper_bound(atom,upper,control)
	# unbounded corner points?
  LIP_CONST = 10
  ## if atom is unbounded, compute a bound based of the paper; take the closest point to the average po
	if !is_bounded(atom)
    if is_fully_unbounded(atom)
      return Inf
    end
    closest = closest_point(atom.A, atom.corner_points)
		return atom.P*(Sprocket.evaluate(upper,closest[1]*control)+closest[2]*LIP_CONST)
	end

  ## if atom is bounded, do a normal edmund madansky upper bound for the problem
	sum = 0.0
	for point in atom.corner_points
		sum+= atom.corner_weights[point]*Sprocket.evaluate(upper,point*control)
	end
	return sum*atom.P
end

function closest_point(point::Sprocket.Point, points::Vector{Sprocket.Point})
  out = collect(zip(points,map(x -> Sprocket.l_one_norm(x - point), points)))

  closest_point = reduce(out) do x, y
    if x[2] < y[2]
      return x
    else
      return y
    end
  end
  return closest_point
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
		out = out && point[i]
	end
	return out
end

function Base.any(point::Sprocket.Point)
	out = false
	for i in eachindex(point)
		out = out || point[i]
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

function add_atom!(c,a::Baucke.Atom) end
function remove_atom!(c,a::Baucke.Atom) end
function get_control!(c) end
function add_cut!(c,cut) end

function reltol_lower(upper::Float64, lower::Float64)
  (upper - lower)/lower
end

function reltol_upper(upper::Float64, lower::Float64)
  (upper - lower)/upper
end

function abstol(upper::Float64, lower::Float64)
  upper - lower
end

end


# TODO
# add way of computing upper bound for unbounded atom
