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

function generate_cut(c_oracle,point::Point)
	(val,grad) = c_oracle(point)
	return Cut(point,val,grad)
end
