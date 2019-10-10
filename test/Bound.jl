using JuMP
using GLPK

include("../src/Variable.jl")
include("../src/Cut.jl")
include("../src/Bound.jl")


hello = Variable(:xi,(),Random())

function my_cutting_plane_oracle()
	function zero_order(vars)
		return vars[:xi]^2
		# return 0
	end
	function first_order(vars)
		grad = deepcopy(vars)
		grad[:xi] = 2*vars[:xi]
		return grad
	end
	return (zero_order,first_order)
end

my_lower = Lowerbound(hello,-10.0)
my_upper = Upperbound(hello, 10.0, 99.9)

origin = Point(Dict(hello => 0.0))
unit = Point(Dict(hello => 1.0))

evaluate(my_lower,unit)

some_cut = generate_cut(my_cutting_plane_oracle(),origin)

update!(my_lower,generate_cut(my_cutting_plane_oracle(),unit))
update!(my_upper,generate_cut(my_cutting_plane_oracle(),unit))

evaluate(my_lower,2*unit)
evaluate(my_upper,2*unit)

update!(my_lower,generate_cut(my_cutting_plane_oracle(),2*unit))
update!(my_upper,generate_cut(my_cutting_plane_oracle(),2*unit))

evaluate(my_lower,2*unit)
evaluate(my_upper,2*unit)

evaluate(my_lower,0*unit)
evaluate(my_upper,0*unit)

update!(my_lower,generate_cut(my_cutting_plane_oracle(),0*unit))
update!(my_upper,generate_cut(my_cutting_plane_oracle(),0*unit))

evaluate(my_lower,0*unit)
evaluate(my_upper,0*unit)

evaluate(my_lower,0.5*unit)
evaluate(my_upper,0.5*unit)
