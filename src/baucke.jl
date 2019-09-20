####
# In this file we define the BauckeAlgorithm which will solve the optimisation problem.
####

# using Sprocket
include("./Sprocket.jl")

function BauckeAlgorithm()
	function initialise(prob)
		lower = Lowerbound(prob,-99)
		upper = Upperbound(prob,99,10)
		return (lower=lower,upper=upper,atoms=atoms,crit=crit)
	end
	function iterate(state)
		atom_with_largest_bound_gap(atoms)
		find_center_point(atom)
		split(atoms,center_point)
		generate_cut(center_point)
		update(lower,cut)
		update(upper,cut)
	end
	function hasmet(crit::Criteria)
	end
	return (intialise,iterate,hasmet)
end

function InitialAtom(prob)
	# get the random structure from the problem
	random = prob[1]()[2]
end
