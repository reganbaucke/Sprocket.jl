####
# In this file, we will implement the functions required for updating and evalutating bounding functions
####
using JuMP
using GLPK #free and open source solver


###
# The two types of bouding functions, we make them distinct because we would like to dispatch different methods on them
###
abstract type Bound end

struct Lowerbound <: Bound
	model::JuMP.Model
	cuts::Array
	queries::Int
end

struct Upperbound <: Bound
	model::JuMP.Model
	cuts::Array
	queries::Int
end


###
# Constructors
###
function Lowerbound(problem, initial_lower)
	# initialise empty model
	model = JuMP.Model(with_optimizer(GLPK.Optimizer))

	## call get the structure of the problem
	structure = problem[1]()

	## for the control variables, set up the JuMP Variables

	for key in keys(structure[1])
		create_jump_variable(model,key,structure[1][key])
	end

	## for the random variables, set up the JuMP Variables
	for key in keys(structure[2])
		create_jump_variable(model,key,structure[2][key])
	end

	## add the epigraph variable
	epi = @variable(model,lower_bound=initial_lower)
	@objective(model,Min, epi)

	cuts = []
	return Lowerbound(model,cuts,0)
end

function Upperbound(problem, initial_lower, lipschitz_bound)
	# initialise empty model
	model = JuMP.Model(with_optimizer(GLPK.Optimizer))

	## call get the structure of the problem
	structure = problem[1]()

	## for the control variables, set up the JuMP Variables
	for key in keys(structure[1])
		create_jump_variable(model,key,structure[1][key])
	end

	## for the random variables, set up the JuMP Variables
	for key in keys(structure[2])
		create_jump_variable(model,key,structure[2][key])
	end

	## add the intercept variable
	intercept = @variable(model,lower_bound=initial_lower)
	@objective(model,Max, intercept)

	cuts=[]

	return Upperbound(model,cuts,0)
end

###
# Evaluate function
###
# point is a tuple of dictionary of symbols to arrays of floats
function evaluate(lower::Lowerbound,point)
	# our goal is to set the variables in the model of lower to the values in point and
	# then solve the model

	control = point[1]
	random = point[2]

	if !isempty(control)
		for key in keys(control)
			if typeof(control[key]) <: Number
				JuMP.fix(lower.model.obj_dict[key],control[key])
			else
				for (index,value) in enumerate(control[key])
					JuMP.fix(lower.model.obj_dict[key][index],control[key][index])
				end
			end
		end
	end

	if !isempty(random)
		for key in keys(random)
			if (typeof(random[key]) <: Number)
				JuMP.fix(lower.model.obj_dict[key],random[key])
			else
				for (index,value) in enumerate(random[key])
					JuMP.fix(lower.model.obj_dict[key][index],random[key][index])
				end
			end
		end
	end

	obj = NaN
	optimize!(lower.model)
	# if termination_status(lower.model) == 1
	# need to perform check to see if model got optimized correctly
	obj = objective_value(lower.model)
	# end

	for key in keys(control)
		if typeof(control[key]) <: Number
			JuMP.unfix(lower.model.obj_dict[key])
		else
			for (index,value) in enumerate(control[key])
				JuMP.unfix(lower.model.obj_dict[key][index])
			end
		end
	end

	for key in keys(random)
		if typeof(random[key]) <: Number
			JuMP.unfix(lower.model.obj_dict[key])
		else
			for (index,value) in enumerate(random[key])
				JuMP.unfix(lower.model.obj_dict[key][index])
			end
		end
	end

	return obj
end

function evaluate(upper::Upperbound,point)
	control = point[1]
	random = point[2]

	# we have to go through and reset the object coefficents for all the variables
	# use JuMP.set_objective_coefficient

	if !isempty(control)
		for key in keys(control)
			if typeof(control[key]) <: Number
				JuMP.set_objective_coefficient(lower.model.obj_dict[key],control[key])
			else
				for (index,value) in enumerate(control[key])
					JuMP.set_objective_coefficient(lower.model,lower.model.obj_dict[key][index],control[key][index])
				end
			end
		end
	end

	if !isempty(random)
		for key in keys(random)
			if (typeof(random[key]) <: Number)
				JuMP.set_objective_coefficient(lower.model.obj_dict[key],random[key])
			else
				for (index,value) in enumerate(random[key])
					JuMP.set_objective_coefficient(lower.model.obj_dict[key][index],random[key][index])
				end
			end
		end
	end

	obj = NaN
	optimize!(lower.model)
	# if termination_status(lower.model) == 1
	# need to perform check to see if model got optimized correctly
	obj = objective_value(lower.model)
	# end

	##
	# Doesn't make sense it ``unset'' the variables in the lower bound case
	##
end

####
# update function
####
function update(lower::Lowerbound,cut)
	# we need to add a constraint to the optimization problem
	cut.obj
	# build up JuMP to use in the constraint
	ex = JuMP.AffExpr(0.0)

	if !isempty(control)
		for key in keys(control)
			if typeof(control[key]) <: Number
				JuMP.set_objective_coefficient(lower.model.obj_dict[key],control[key])
			else
				for (index,value) in enumerate(control[key])
					JuMP.set_objective_coefficient(lower.model,lower.model.obj_dict[key][index],control[key][index])
				end
			end
		end
	end

	if !isempty(random)
		for key in keys(random)
			if (typeof(random[key]) <: Number)
				JuMP.set_objective_coefficient(lower.model.obj_dict[key],random[key])
			else
				for (index,value) in enumerate(random[key])
					JuMP.set_objective_coefficient(lower.model.obj_dict[key][index],random[key][index])
				end
			end
		end
	end


	@constraint(lower.model,cut.obj)
end

function update(upper::Upperbound,cut)

end



# JuMP hacks
# i really shouldn't have to write these myself
function create_jump_variable(model::JuMP.Model,name::Symbol,size)
	# if a scalar
	if size == ()
		# my_var = JuMP.variable_type(model)(JuMP.undef)
		my_var = JuMP.add_variable(model, JuMP.build_variable(x-> (), JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)), JuMP.string(String(name)))
		(JuMP.object_dictionary(model))[name] = my_var
		return my_var
	end

	# my_var = JuMP.Array{JuMP.variable_type(model)}(JuMP.undef, (JuMP.length(Base.OneTo(9)),)...)
	my_var = JuMP.Array{JuMP.variable_type(model)}(JuMP.undef, size...)
	for i in CartesianIndices(my_var)
		my_var[i] = JuMP.add_variable(model, JuMP.build_variable(x-> (), JuMP.VariableInfo(false, NaN, false, NaN, false, NaN, false, NaN, false, false)), JuMP.string(String(name), "[", string(Tuple(i))  , "]"))
	end
	(JuMP.object_dictionary(model))[name] = my_var
	return my_var
end

function do_for_each(point)
end
