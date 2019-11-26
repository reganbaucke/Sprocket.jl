struct Cut
	point::Point
	value::Number
	grad::Point
end

# function generate_cut(c_oracle,point::Point)
# 	(zero,first) = c_oracle
# 	val = zero(point)
# 	grad = first(point)
# 	return Cut(point,val,grad)
# end

function generate_cut(c_oracle, point::Point)
	(val,grad) = c_oracle(point)
	return Cut(point,val,grad)
end

function apply_cut_at(cut::Cut, point::Point)
	new_point = sub_point(cut.point,setdiff(get_vars(cut.point),get_vars(point)))
	new_value = cut.value + inner(point-sub_point(cut.point,get_vars(point)), sub_point(cut.grad,get_vars(point)))
	new_grad = sub_point(cut.grad,setdiff(get_vars(cut.grad),get_vars(point)))
	Cut(new_point,new_value,new_grad)
end



###
# Small suite of tests
###
function test()

end
