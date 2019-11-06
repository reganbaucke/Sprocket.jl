####
# In this file we define the Sample Average Approximation Algorithm (SAA) which will solve the optimisation problem exactly
####
module SAA

using ..Sprocket
using JuMP
using GLPK

abstract type AlgorithmVersion end
struct MultiCut <: AlgorithmVersion  end
struct AverageCut <: AlgorithmVersion end

struct State
   lower_bound
   criteria::Sprocket.Criteria
end

function SAAAlgorithm(samples, version::AlgorithmVersion)
   function initialise(prob)
      ###
      # Set up the control problem
      ###
      lower_bound = build_lower_problem(prob,samples,version)

      ##
      # Set up the initial status of the criteria
      ###
      new_criteria = Sprocket.initial_criteria()

      return SAA.State(lower_bound,new_criteria)
   end

   function iterate(state,prob)
      (control, lower_bound_val) = get_new_control(state.lower_bound)
      # compute cuts from oracle at points

      if version == AverageCut()
         # average the cuts
      elseif version == MultiCut()
         for sample in samples
            cut = Sprocket.generate_cut(prob.c_oracle,sample*control)
            expression = gen_constraint_from_cut(state.lower_bound, sample, cut)
            apply_constraint(state.lower_bound, sample, expression)
            add_cut!(state.lower_bound, cut)
         end
      end

      # update criteria
      new_criteria = Sprocket.update_with(state.criteria,
         iterations = state.criteria.iterations + 1,
         lower = lower_bound_val
      )

      return SAA.State(state.lower_bound, new_criteria)
   end

   function hasmet(crit::Sprocket.Criteria,state)
      Sprocket.met(state.criteria,crit)
   end

   return (initialise,iterate,hasmet)
end

struct ControlProblem
   model::JuMP.Model
   control_vars::Sprocket.Point
   epi_vars::Dict{Sprocket.Point,JuMP.VariableRef}
   cut_constraints::Dict{Sprocket.Point,Vector{JuMP.ConstraintRef}}
   cuts::Vector{Sprocket.Cut}
end

function build_lower_problem(prob::Sprocket.Problem, samples::Vector{Sprocket.Point}, version::SAA.MultiCut)
   # determine the control variables
   controls = filter(x->x.type == Sprocket.Control(),prob.vars)
   if isempty(controls)
      return nothing
   end

   my_model = copy(prob.model)
   @objective(my_model,Min,0)

   epi_vars = Dict{Sprocket.Point,JuMP.VariableRef}()
   for sample in samples
      epi_vars[sample] = @variable(my_model, lower_bound = -99)
      JuMP.set_objective_coefficient(my_model,epi_vars[sample],1.0/length(samples))
   end

   cut_constraints = Dict{Sprocket.Point, Array{JuMP.ConstraintRef}}()
   for sample in samples
      cut_constraints[sample] = JuMP.ConstraintRef[]
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

function gen_constraint_from_cut(prob, sample, cut)
   reduced_cut = Sprocket.apply_cut_at(cut, sample)

   # build up constraint in an affine expression
   ex = JuMP.AffExpr(0.0)

   for i in eachindex(prob.control_vars)
      JuMP.add_to_expression!(ex,reduced_cut.grad[i]*(prob.control_vars[i] - reduced_cut.point[i]))
   end

   # add the constant term of the cut value
   JuMP.add_to_expression!(ex,reduced_cut.value)

   # add the epigraph variable associated with this atom
   JuMP.add_to_expression!(ex, - prob.epi_vars[sample])
   return ex
end

function apply_constraint(prob::ControlProblem, sample, ex)
   ref = @constraint(prob.model, 0 >= ex)
   push!(prob.cut_constraints[sample],ref)
end

function add_cut!(prob::ControlProblem,cut::Sprocket.Cut)
   push!(prob.cuts,cut)
end


function test()
   function build_problem()
      #create an empty problem
      problem = Sprocket.Problem()

      # add one real valued random variable to the problem
      xi = Sprocket.Variable(name=:xi, size=(),type=Sprocket.Random())
      Sprocket.add_variable(problem,xi)

      #set the domain of the random variable
      Sprocket.set_domain(problem,(Sprocket.Point(Dict(xi =>0.0)), Sprocket.Point(Dict(xi => 1.0))))


      u = Sprocket.Variable(name=:u, size=(),type=Sprocket.Control())
      Sprocket.add_variable(problem,u)
      @variable(problem.model, -1.0 <= u <= 1.0)

      # objective function is defined through the cutting plane oracle
      function my_cutting_plane_oracle()
         function oracle(vars)
            value = vars[:xi]^2 + 0.1*vars[:u]*vars[:xi] + vars[:u]^2
            grad = deepcopy(vars)
            grad[:xi] = 2.0*vars[:xi] + 0.1*vars[:u]
            grad[:u]  = 0.1*vars[:xi] + 2.0*vars[:u]
            return (value,grad)
         end
         return oracle
      end

      # distribution of random variables is defined through the measure oracle
      function my_measure_oracle()
         function oracle(vars_1,vars_2)
            left = vars_1[:xi]
            right = vars_2[:xi]
            if vars_1[:xi] <= 0
               left = 0.0
            end
            if vars_1[:xi] >= 1
               left = 1.0
            end
            if vars_2[:xi] <= 0
               right = 0.0
            end
            if vars_2[:xi] >= 1
               right = 1.0
            end
            fresh = deepcopy(vars_1)
            fresh[:xi] = (right^2 - left^2)/2
            return (right-left,fresh)
         end
         return oracle
      end

      problem.c_oracle = my_cutting_plane_oracle()
      problem.m_oracle = my_measure_oracle()

      return problem
   end

   my_prob = build_problem()
   this_point = Sprocket.Point(Dict(my_prob.vars[:xi] => 1.0, my_prob.vars[:u] => 4.0))
   random_part = Sprocket.sub_point(this_point,my_prob.vars[:xi])
   before = Sprocket.generate_cut(my_prob.c_oracle,this_point)
   hello = Sprocket.generate_cut(my_prob.c_oracle,this_point) |> x -> Sprocket.apply_cut_at(x,0*Sprocket.sub_point(this_point,my_prob.vars[:xi]))


   samples = Sprocket.Point[]
   for i = 1 : 10
      push!(samples, Sprocket.random_sample(my_prob.m_oracle,10.0, my_prob.domain))
   end

   my_control_problem = build_control_problem(my_prob, samples, SAA.MultiCut())
end

end
