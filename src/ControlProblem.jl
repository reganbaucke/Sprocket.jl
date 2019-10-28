# include("./Sprocket.jl")
# include("./Baucke.jl")

using JuMP
using GLPK
using .Sprocket
using .Baucke

struct ControlProblem
	model::JuMP.Model
	control_vars::Sprocket.Point
	epi_vars::Dict{Baucke.Atom,JuMP.VariableRef}
	cut_constraints::Dict{Baucke.Atom,Array{JuMP.ConstraintRef}}
	cuts::Vector{Sprocket.Cut}
end

function build_control_problem(prob::Sprocket.Problem)
	# determine the control variables
	controls = filter(x->x.type == Sprocket.Control(),prob.vars)
	if isempty(controls)
		return nothing
	end

	my_model = copy(prob.model)
	@objective(my_model,Min,0)

	epi_vars = Dict{Baucke.Atom,JuMP.VariableRef}()
	cut_constraints = Dict{Baucke.Atom,Array{JuMP.ConstraintRef}}()


	storage = Dict{Sprocket.Variable,Any}()
	for control in controls
		storage[control] = my_model.obj_dict[control.name]
	end

	cuts = Sprocket.Cut[]
	Some(ControlProblem(my_model,Sprocket.Point(storage),epi_vars,cut_constraints,cuts))
end

function add_atom(prob::ControlProblem,atom::Baucke.Atom)
	# create a new jump variable
	anon = @variable(prob.model,lower_bound = -99)
	prob.epi_vars[atom] = anon
	prob.cut_constraints[atom] = JuMP.ConstraintRef[]

	# set objective value of the as the probability of the atom
	JuMP.set_objective_coefficient(prob.model, anon,atom.P)
end

function gen_constraints_for_new_atom(prob::ControlProblem,atom::Baucke.Atom)
	map(x -> gen_constraint_from_cut(prob,atom,x),prob.cuts) |> (x -> map(y -> apply_constraint(prob,atom,y),x))
end

function update_all_atoms_with_cut(prob::ControlProblem,cut)
	for atom in keys(prob.epi_vars)
		gen_constraint_from_cut(prob,atom,cut) |> x-> apply_constraint(prob,atom,x)
	end
end

function delete_atom(prob::ControlProblem,atom::Baucke.Atom)
	# delete constraints associated to that atom from prob.model
	map(x -> JuMP.delete(prob.model,x),prob.cut_constraints[atom])

	# delete variable associated to that atom from prob.model
	JuMP.delete(prob.model,prob.epi_vars[atom])

	# delete atom from prob.epi_vars Dict
	delete!(prob.epi_vars,atom)
end

function gen_constraint_from_cut(prob::ControlProblem, atom::Baucke.Atom, cut::Sprocket.Cut)
	# first reduce the cut the average position
	reduced_cut = Sprocket.apply_cut_at(cut,atom.A)

	# build up JuMP to use in the constraint
	ex = JuMP.AffExpr(0.0)

	# @assert(get_vars(reduced_cut.point) == get_vars(vars))

	for i in eachindex(prob.control_vars)
		JuMP.add_to_expression!(ex,reduced_cut.grad[i]*(prob.control_vars[i] - reduced_cut.point[i]))
	end

	# add the constant term of the cut value
	JuMP.add_to_expression!(ex,reduced_cut.value)

	# add the epigraph variable associated with this atom
	JuMP.add_to_expression!(ex, - prob.epi_vars[atom])

	return ex
end

function get_new_control(prob::ControlProblem)
	if prob.model.moi_backend.optimizer == nothing
		JuMP.optimize!(prob.model,with_optimizer(GLPK.Optimizer))
	else
		JuMP.optimize!(prob.model)
	end
	map(JuMP.value,prob.control_vars)
end

function apply_constraint(prob::ControlProblem,atom::Baucke.Atom,ex)
	ref = @constraint(prob.model, 0 >= ex)
	push!(prob.cut_constraints[atom],ref)
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
		Sprocket.set_domain(problem,Baucke.rect_hull,(Sprocket.Point(Dict(xi =>0.0)), Sprocket.Point(Dict(xi => 1.0))))


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
	# println(before)
	println(before.value)
	println(before.point)
	println(before.grad)
	println("-----------")
	println(0*Sprocket.sub_point(this_point,my_prob.vars[:xi]))
	println("-----------")
	# println(hello)
	println(hello.value)
	println(hello.point)
	println(hello.grad)


	my_control_problem = build_control_problem(my_prob)

	atom = Baucke.Atom()
	atom.corner_points = my_prob.domain

	atom.P = Baucke.compute_probability(atom,my_prob.m_oracle)
	atom.A = Baucke.compute_average_point(atom,my_prob.m_oracle)

	Baucke.compute_weights!(atom)

	add_atom(my_control_problem,atom)
	println(my_control_problem)

	println(atom.A)
	point_to_sample =  atom.A*Sprocket.Point(Dict(my_prob.vars[:u] => 0.0))
	# println(Sprocket.generate_cut(my_prob.c_oracle,atom.A*Sprocket.Point(Dict(my_prob.vars[:u] => 0.0))))

	# GenerateConstraintFromCut(my_control_problem,atom,Sprocket.generate_cut(my_prob.c_oracle,atom.A))
	this = gen_constraint_from_cut(my_control_problem,atom,Sprocket.generate_cut(my_prob.c_oracle,point_to_sample))
	println(this)
	apply_constraint(my_control_problem,atom,this)
end
