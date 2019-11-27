####
# In this file we define the ExactAlgorithm which will solve the optimisation problem exactly
####
module KuhnMaggioni

using Combinatorics
using ..Sprocket
using JuMP
using GLPK

#structs for use
mutable struct Atom
   corner_points
   corner_weights
   P
   A
   function Atom()
      return new()
   end
end

# include("./ControlProblem.jl")

struct State
   lower
   atoms::Set{KuhnMaggioni.Atom}
   criteria::Sprocket.Criteria
end

function Algorithm(divisions)
   function initialise(prob)
      # verify that the domain is bounded

      # generate lists of pairs of points to generate atoms
      pairs = cartesian_product(divisions,prob.domain)

      # divide up the domain and create all the atoms according to the division parameter
      atoms = Set()

      for pair in pairs
         atom = KuhnMaggioni.Atom()

         atom.corner_points = rect_hull(pair...)
         atom.P = compute_probability(atom,prob.m_oracle)
         atom.A = compute_average_point(atom,prob.m_oracle)
         compute_weights!(atom)

         push!(atoms,atom)
      end

      lower_bound = build_lower_problem(prob,atoms)

      # start the criteria with the number of measure oracle calls made
      criteria = Sprocket.update_with(Sprocket.initial_criteria(),
         m_calls = Sprocket.fold( *, 1, divisions)
      )

      return KuhnMaggioni.State(lower_bound, atoms, criteria)
   end


   function iterate(state,prob)
      (control, lower_bound_val) = get_new_control(state.lower)

      ###
      # Update compute an upper bound
      ###
      upper_bound_val = compute_upper_bound(prob.c_oracle,state.atoms,control)

      # compute cuts from oracle at points
      for atom in state.atoms
         cut = Sprocket.generate_cut(prob.c_oracle,atom.A*control)
         expression = gen_constraint_from_cut(state.lower, atom, cut)
         apply_constraint(state.lower, atom, expression)
         add_cut!(state.lower, cut)
      end

      # update criteria
      new_criteria = Sprocket.update_with(state.criteria,
         iterations = state.criteria.iterations + 1,
         lower = lower_bound_val,
         upper = upper_bound_val,
         abstol = abstol(upper_bound_val,lower_bound_val),
         reltol = reltol_upper(upper_bound_val,lower_bound_val)
      )

      # compute the bound gap at the new control
      # biggest_bound = largest_bound_gap(state.atoms,(state.lower,state.upper))

      # return new_state
      return KuhnMaggioni.State(state.lower, state.atoms, new_criteria)
   end

   function hasmet(crit::Sprocket.Criteria,state)
      Sprocket.met(state.criteria,crit)
   end

   return (initialise,iterate,hasmet)
end


function compute_probability(atom::KuhnMaggioni.Atom, oracle)
   return oracle(get_generating_pair(atom)...)[1]
end

function compute_average_point(atom::KuhnMaggioni.Atom, oracle)
   return oracle(get_generating_pair(atom)...)[2]/oracle(get_generating_pair(atom)...)[1]
end

function cartesian_product(division,domain)
   #
   # assert that the domain is bounded
   #

   # instead of two points, make one list of pairs of points
   new_domain = zip(domain...)
   # println(new_domain)


   # make function for use in combine
   function merge(x,y)
      one_d_atoms = []
      for j in 1:x
         push!(one_d_atoms,
            (y[1] + (j - 1)*(y[2] - y[1])/x,
            y[1] + (j)*(y[2] - y[1])/x))
      end
      one_d_atoms
   end

   point_of_pairs = Sprocket.combine(merge,division,new_domain)

   # we now have a point of lists of pairs
   # we want to have a list of pairs of points
   list = []
   for var in eachindex(point_of_pairs)
      inner_list = []
      for pair in point_of_pairs[var]
         push!(inner_list, (Sprocket.Point(Dict( var => pair[1])), Sprocket.Point(Dict( var => pair[2]))) )
      end
      push!(list, inner_list)
   end

   function cartesian(listOfPoints,newList)
      out = []
      if isempty(listOfPoints)
         return newList
      end
      for orig in listOfPoints
         for new in newList
            push!(out, (orig[1]*new[1], orig[2]*new[2]))
            # push!(out, orig[2]*new[1])
            # push!(out, orig[1]*new[2])
            # push!(out, orig[2]*new[2])
         end
      end
      out
   end

   out = []
   for inner_list in list
      out = cartesian(out, inner_list)
   end
   out

end

struct ControlProblem
   model::JuMP.Model
   control_vars::Sprocket.Point
   epi_vars::Dict{KuhnMaggioni.Atom,JuMP.VariableRef}
   cut_constraints::Dict{KuhnMaggioni.Atom,Vector{JuMP.ConstraintRef}}
   cuts::Vector{Sprocket.Cut}
end

function build_lower_problem(prob::Sprocket.Problem, atoms)
   # determine the control variables
   controls = filter(x->x.type == Sprocket.Control(),prob.vars)

   my_model = copy(prob.model)
   @objective(my_model,Min,0)

   epi_vars = Dict{KuhnMaggioni.Atom,JuMP.VariableRef}()
   for atom in atoms
      epi_vars[atom] = @variable(my_model, lower_bound = -99)
      JuMP.set_objective_coefficient(my_model,epi_vars[atom],atom.P)
   end

   cut_constraints = Dict{KuhnMaggioni.Atom, Array{JuMP.ConstraintRef}}()
   for atom in atoms
      cut_constraints[atom] = JuMP.ConstraintRef[]
   end

   storage = Dict{Sprocket.Variable,Any}()
   for control in controls
      storage[control] = my_model.obj_dict[control.name]
   end

   cuts = Sprocket.Cut[]

   ControlProblem(my_model,Sprocket.Point(storage),epi_vars,cut_constraints,cuts)
end

function get_new_control(prob::ControlProblem)
   if prob.model.moi_backend.optimizer == nothing
      JuMP.optimize!(prob.model,with_optimizer(GLPK.Optimizer))
   else
      JuMP.optimize!(prob.model)
   end
   (map(JuMP.value,prob.control_vars), JuMP.objective_value(prob.model))
end

function compute_upper_bound(c_oracle,atoms,control)
   oracle_cache = Dict()
   upper = 0.0
   for atom in atoms
      upper_on_this_atom = 0.0
      for point in atom.corner_points
         if point in keys(oracle_cache)
            upper_on_this_atom += oracle_cache[point]*atom.corner_weights[point]*atom.P
         else
            oracle_cache[point] = c_oracle(point*control)[1]
            upper_on_this_atom += oracle_cache[point]*atom.corner_weights[point]*atom.P
         end
      end
      upper += upper_on_this_atom
   end

   upper
end

function compute_weights!(atom::KuhnMaggioni.Atom)
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

function lower_bound(atom,lower,control)
   atom.P*(Sprocket.evaluate(lower,atom.A*control))
end

function is_bounded(atom::KuhnMaggioni.Atom)
   for point in atom.corner_points
      if any(map(x -> x == Inf || x == -Inf, point))
         return false
      end
   end
   return true
end

function is_fully_unbounded(atom::KuhnMaggioni.Atom)
   for point in get_generating_pair(atom)
      if !all(map(x -> x == Inf || x == -Inf, point))
         return false
      end
   end
   return true
end

function gen_constraint_from_cut(prob, atom, cut)
   reduced_cut = Sprocket.apply_cut_at(cut, atom.A)

   # build up constraint in an affine expression
   ex = JuMP.AffExpr(0.0)

   for i in eachindex(prob.control_vars)
      JuMP.add_to_expression!(ex,reduced_cut.grad[i]*(prob.control_vars[i] - reduced_cut.point[i]))
   end

   # add the constant term of the cut value
   JuMP.add_to_expression!(ex,reduced_cut.value)

   # add the epigraph variable associated with this atom
   JuMP.add_to_expression!(ex, - prob.epi_vars[atom])
   return ex
end

function apply_constraint(prob::ControlProblem, atom, ex)
   ref = @constraint(prob.model, 0 >= ex)
   push!(prob.cut_constraints[atom],ref)
end

function add_cut!(prob::ControlProblem,cut::Sprocket.Cut)
   push!(prob.cuts,cut)
end

function get_probability(m_oracle,a::KuhnMaggioni.Atom)
   (zero,first) = m_oracle
   return zero(get_generating_pair(a...))/zero(get_generating_pair(a...))
end

function get_probability(m_oracle,a::KuhnMaggioni.Atom)
   (zero,first) = m_oracle
   return zero(get_generating_pair(a...))/zero(get_generating_pair(a...))
end

function get_average_point(m_oracle,a::KuhnMaggioni.Atom)
   (zero,first) = m_oracle
   return first(get_generating_pair(a...))
end

function get_generating_pair(a::KuhnMaggioni.Atom)
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


function Base.string(atom::KuhnMaggioni.Atom)
   out = ""
   for point in atom.corner_points
      out *= string(point) * "\n"
   end

   out*= "probability: $(atom.P) \n"
   out*= "barycenter: $(atom.A)"
end

function Base.show(io::IO,atom::KuhnMaggioni.Atom)
   out = ""

   for point in atom.corner_points
      out *= string(point) * "\n"
   end

   out*= "probability: $(atom.P) \n"
   out*= "barycenter: $(atom.A)"
   println(io,out)
end

function reltol_lower(upper::Float64, lower::Float64)
   (upper - lower)/lower
end

function reltol_upper(upper::Float64, lower::Float64)
   (upper - lower)/upper
end

function abstol(upper::Float64, lower::Float64)
   upper - lower
end

function list_of_cartesian_product(divisions::Dict)
   output = []
   for keys in divisions

   end
end

function reltol_lower(upper::Float64, lower::Float64)
   abs(upper - lower)/abs(lower)
end

function reltol_upper(upper::Float64, lower::Float64)
   abs(upper - lower)/abs(upper)
end

function abstol(upper::Float64, lower::Float64)
   upper - lower
end

end
