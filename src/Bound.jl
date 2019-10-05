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
	vars
	epi::JuMP.VariableRef
	model::JuMP.Model
	cuts::Array
	queries::Int
end

struct Upperbound <: Bound
	vars
	model::JuMP.Model
	cuts::Array
	queries::Int
end


###
# Constructors
###
function Lowerbound(vars, initial_lower)
	# initialise empty model
	model = JuMP.Model(with_optimizer(GLPK.Optimizer))


	## call get the structure of the problem

	## for the control variables, set up the JuMP Variables
	for key in keys(vars)
		create_jump_variable(model,key,vars[key][2])
	end

	## add the epigraph variable
	epi = @variable(model,lower_bound=initial_lower)
	@objective(model,Min, epi)

	cuts = []
	return Lowerbound(copy(vars),epi,model,cuts,0)
end

function Upperbound(vars, initial_upper, lipschitz_bound)
	# initialise empty model
	model = JuMP.Model(with_optimizer(GLPK.Optimizer))


	## for the control variables, set up the JuMP Variables
	for key in keys(vars)
		create_jump_variable(model,key,vars[key][2],(-lipschitz_bound,lipschitz_bound))
	end

	## add the intercept variable
	@variable(model,intercept,upper_bound=initial_upper)
	@objective(model,Max, intercept)

	cuts=[]

	return Upperbound(copy(vars),model,cuts,0)
end

###
# Evaluate function
###
function evaluate(lower::Lowerbound,vars)
	# our goal is to set the variables in the model of lower to the values in point and
	# then solve the model

	# check that we are sending in the right points
	@assert(keys(vars)==keys(lower.vars))


	for key in keys(vars)
		if typeof(vars[key]) <: Number
			JuMP.fix(lower.model.obj_dict[key],vars[key])
		else
			for (index,value) in enumerate(vars[key])
				JuMP.fix(lower.model.obj_dict[key][index],value)
			end
		end
	end

	obj = NaN
	optimize!(lower.model)
	# if termination_status(lower.model) == 1
	# need to perform check to see if model got optimized correctly
	obj = objective_value(lower.model)
	# end

	for key in keys(vars)
		if typeof(vars[key]) <: Number
			JuMP.unfix(lower.model.obj_dict[key])
		else
			for (index,value) in enumerate(vars[key])
				JuMP.unfix(lower.model.obj_dict[key][index])
			end
		end
	end

	return obj
end

function evaluate(upper::Upperbound,vars)
	# for code readability
	@assert(keys(vars)==keys(upper.vars))

	# we have to go through and reset the object coefficents for all the variables
	# use JuMP.set_objective_coefficient

	for key in keys(vars)
		if typeof(vars[key]) <: Number
			JuMP.set_objective_coefficient(upper.model,upper.model.obj_dict[key],vars[key])
		else
			for (index,value) in enumerate(vars[key])
				JuMP.set_objective_coefficient(upper.model,upper.model.obj_dict[key][index],vars[key][index])
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

	@assert(keys(cut.point) == keys(lower.vars))

	for key in keys(lower.vars)
		if typeof(cut.point[key]) <: Number
			JuMP.add_to_expression!(ex,cut.grad[key]*(lower.model.obj_dict[key] - cut.point[key]))
		else
			for (index,value) in enumerate(cut.point[key])
				JuMP.add_to_expression!(ex,cut.grad[key][index]*(lower.model.obj_dict[key][index] - cut.point[key][index]))
			end
		end
	end

	return @constraint(lower.model,lower.:epi >= cut.value +  ex)
end

function update!(upper::Upperbound,cut)
	# push the cut to the list of cuts
	push!(upper.cuts,cut)

	@assert(keys(cut.point)==keys(upper.vars))

	# build up JuMP to use in the constraint
	ex = JuMP.AffExpr(0.0)

	for key in keys(upper.vars)
		if typeof(cut.point[key]) <: Number
			JuMP.add_to_expression!(ex,cut.point[key],upper.model.obj_dict[key])
		else
			for (index,value) in enumerate(cut.point[key])
				JuMP.add_to_expression!(ex,cut.point[key][index],upper.model.obj_dict[key][index])
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


# TODO
# make a new entry in the UPPERBOUND struct that keeps track of the intercept variable, like
# what is done for the `epi' variable in the LOWERBOUND.
