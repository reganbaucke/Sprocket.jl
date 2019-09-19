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
	for key in keys(structure.control)
		create_jump_variable(model,key,structure.control[key])
	end

	## for the random variables, set up the JuMP Variables
	for key in keys(structure.random)
		create_jump_variable(model,key,structure.random[key])
	end

	## add the epigraph variable
	@variable(model,epi,lower_bound=initial_lower)
	@objective(model,Min, epi)

	cuts = []
	return Lowerbound(model,cuts,0)
end

function Upperbound(problem, initial_upper, lipschitz_bound)
	# initialise empty model
	model = JuMP.Model(with_optimizer(GLPK.Optimizer))

	## call get the structure of the problem
	structure = problem[1]()

	## for the control variables, set up the JuMP Variables
	for key in keys(structure.control)
		create_jump_variable(model,key,structure.control[key],(-lipschitz_bound,lipschitz_bound))
	end

	## for the random variables, set up the JuMP Variables
	for key in keys(structure.random)
		create_jump_variable(model,key,structure.random[key],(-lipschitz_bound,lipschitz_bound))
	end

	## add the intercept variable
	@variable(model,intercept,upper_bound=initial_upper)
	@objective(model,Max, intercept)

	cuts=[]

	return Upperbound(model,cuts,0)
end

###
# Evaluate function
###
# point is a tuple of dictionary of symbols to arrays of floats
function evaluate(lower::Lowerbound,point::Point)
	# our goal is to set the variables in the model of lower to the values in point and
	# then solve the model

	# for code readability
	control = point.control
	random = point.random

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
			if !(typeof(random[key]) <: AbstractArray)
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

function evaluate(upper::Upperbound,point::Point)
	# for code readability
	control = point.control
	random = point.random

	# we have to go through and reset the object coefficents for all the variables
	# use JuMP.set_objective_coefficient

	if !isempty(control)
		for key in keys(control)
			if typeof(control[key]) <: Number
				JuMP.set_objective_coefficient(upper.model,upper.model.obj_dict[key],control[key])
			else
				for (index,value) in enumerate(control[key])
					JuMP.set_objective_coefficient(upper.model,lower.model.obj_dict[key][index],control[key][index])
				end
			end
		end
	end

	if !isempty(random)
		for key in keys(random)
			if (typeof(random[key]) <: Number)
				JuMP.set_objective_coefficient(upper.model,upper.model.obj_dict[key],random[key])
			else
				for (index,value) in enumerate(random[key])
					JuMP.set_objective_coefficient(upper.model,upper.model.obj_dict[key][index],random[key][index])
				end
			end
		end
	end

	obj = NaN
	optimize!(upper.model)
	# if termination_status(lower.model) == 1
	# need to perform check to see if model got optimized correctly
	obj = objective_value(upper.model)
	# end

	return obj

	##
	# Doesn't make sense it ``unset'' the variables in the upper bound case
	##
end

####
# update functions
####
function update!(lower::Lowerbound,cut)
	# push the cut to the list of cuts
	push!(lower.cuts,cut)

	# build up JuMP to use in the constraint
	ex = JuMP.AffExpr(0.0)

	if !isempty(cut.grad.control)
		for key in keys(cut.grad.control)
			if typeof(cut.grad.control[key]) <: Number
				JuMP.add_to_expression!(ex,cut.grad.control[key],(lower.model.obj_dict[key] - cut.point.control[key]))
			else
				for (index,value) in enumerate(cut.grad.control[key])
					JuMP.add_to_expression!(ex,cut.grad.control[key][index],(lower.model.obj_dict[key][index] - cut.point.control[key][index]))
				end
			end
		end
	end

	if !isempty(cut.grad.random)
		for key in keys(cut.grad.random)
			if typeof(cut.grad.random[key]) <: Number
				JuMP.add_to_expression!(ex,cut.grad.random[key],(lower.model.obj_dict[key] - cut.point.random[key]))
			else
				for (index,value) in enumerate(cut.grad.random[key])
					JuMP.add_to_expression!(ex,cut.grad.random[key][index],(lower.model.obj_dict[key][index] - cut.point.random[key][index]))
				end
			end
		end
	end

	return @constraint(lower.model,lower.model.obj_dict[:epi] >= cut.value +  ex)
end

function update!(upper::Upperbound,cut)
	# push the cut to the list of cuts
	push!(upper.cuts,cut)

	# build up JuMP to use in the constraint
	ex = JuMP.AffExpr(0.0)

	if !isempty(cut.grad.control)
		for key in keys(cut.grad.control)
			if typeof(cut.grad.control[key]) <: Number
				JuMP.add_to_expression!(ex,cut.grad.point[key],upper.obj_dict[key])
			else
				for (index,value) in enumerate(cut.grad.control[key])
					JuMP.add_to_expression!(ex,cut.grad.control[key][index],upper.obj_dict[key][index])
				end
			end
		end
	end

	if !isempty(cut.grad.random)
		for key in keys(cut.grad.random)
			if typeof(cut.grad.random[key]) <: Number
				JuMP.add_to_expression!(ex,cut.point.random[key],upper.model.obj_dict[key])
			else
				for (index,value) in enumerate(cut.point.random[key])
					JuMP.add_to_expression!(ex,cut.point.random[key][index],upper.model.obj_dict[key][index])
				end
			end
		end
	end

	return @constraint(upper.model,upper.model.obj_dict[:intercept] + ex  <= cut.value)
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

function create_jump_variable(model::JuMP.Model,name::Symbol,size,bounds)
	# if a scalar
	if size == ()
		# my_var = JuMP.variable_type(model)(JuMP.undef)
		my_var = JuMP.add_variable(model, JuMP.build_variable(x-> (), JuMP.VariableInfo(true, bounds[1], true, bounds[2], false, NaN, false, NaN, false, false)), JuMP.string(String(name)))
		(JuMP.object_dictionary(model))[name] = my_var
		return my_var
	end

	# my_var = JuMP.Array{JuMP.variable_type(model)}(JuMP.undef, (JuMP.length(Base.OneTo(9)),)...)
	my_var = JuMP.Array{JuMP.variable_type(model)}(JuMP.undef, size...)
	for i in CartesianIndices(my_var)
		my_var[i] = JuMP.add_variable(model, JuMP.build_variable(x-> (), JuMP.VariableInfo(true, bounds[1], true, bounds[2], false, NaN, false, NaN, false, false)), JuMP.string(String(name), "[", string(Tuple(i))  , "]"))
	end
	(JuMP.object_dictionary(model))[name] = my_var
	return my_var
end
