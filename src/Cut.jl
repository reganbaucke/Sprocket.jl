struct Cut
	point
	value::Number
	grad
end

function generate_cut(c_oracle,vars)
	(zero,first) = c_oracle
	val = zero(vars)
	grad = first(vars)
	return Cut(vars,val,grad)
end
