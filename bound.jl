####
# In this file, we will implement the functions required for updating and evalutating bounding functions
####

###
# The two types of bouding functions
###
struct Lowerbound <: Bound end
struct Upperbound <: Bound end

using JuMP
using GLPK

###
# Constructors
###
function Lowerbound(problem) end
function Upperbound(problem) end

###
# Evaluate function
###
function evaluate(lower::Lowerbound,point) end
function evaluate(upper::Upperbound,point) end

####
# update function
####
function update(lower::Lowerbound,cut) end
function update(lower::Lowerbound,cut) end
