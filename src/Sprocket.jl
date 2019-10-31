module Sprocket

include("./utils.jl")
include("./Variable.jl")
include("./Cut.jl")
include("./Problem.jl")
include("./Bound.jl")
include("./Criteria.jl")

function Sprocket.solve(algorithm,prob,crit)
  (init,iterate,check) = algorithm
  states = []

  push!(states,init(prob))
  while !check(crit,states[end])
    push!(states,iterate(states[end],prob))
  end
  return states
end


###
# Generate a random sample from the measure oracle using the acceptance rejection algorithm
###
function random_sample(m_oracle, big_M::Real, domain)
  candidate_point = unit(get_vars(domain[1]))

  while true
    net_density = 1.0

    for i in eachindex(candidate_point)
      (candidate_point[i], density) = generate(domain[1][i], domain[2][i])
      net_density = net_density * density
    end

    uniform_zero_one = rand()

    ### acceptance rejection test
    if uniform_zero_one*big_M*net_density < estimate_density(m_oracle,candidate_point)
      ### acceptance
      return candidate_point
    end
  end
end

###
# estimate the value of the pdf encoded in the measure oracle at a given point.
###
function estimate_density(m_oracle, point)
  step_size = 0.001
  # do central differences on the cdf (given by m_oracle) to get estimate of density
  m_oracle(point - step_size*unit(get_vars(point)), point + step_size*unit(get_vars(point)))[1]/(2*step_size)^(dimension(point))
end

###
# generates a random variate in the between left and right. If left and right are bounded, then the variate is uniform,
# if the left and right are unbounded, generate a standard cauchy variate. We generate a cauchy variate because it has the fatest tails and
# will always dominate any other distribution which has first moments
###
function generate(left::Float64, right::Float64)
  uniform_zero_one = rand()
  if left == -Inf
    if right == Inf
      return (s_cauchy_variate(uniform_zero_one), s_cauchy_variate(uniform_zero_one) |> s_cauchy_pdf)
    else
      return (semi_cauchy_variate(uniform_zero_one,left,right), semi_cauchy_pdf(semi_cauchy_variate(uniform_zero_one,left,right),left,right))
    end
  else
    if right == Inf
      return (semi_cauchy_variate(uniform_zero_one,left,right), semi_cauchy_pdf(semi_cauchy_variate(uniform_zero_one,left,right),left,right))
    else
      return (uniform_variate(uniform_zero_one,left,right), uniform_pdf(uniform_variate(uniform_zero_one,left,right),left,right))
    end
  end
end

function s_cauchy_variate(u::Float64)
  tan(pi*(u - 0.5))
end

function s_cauchy_pdf(x)
  (1/pi)*(1/(1 + x^2))
end

function uniform_variate(u::Float64, left::Float64, right::Float64)
  u*(right - left) + left
end

function uniform_pdf(x::Float64, left::Float64, right::Float64)
  if (left <= x)  && (x <= right)
    return 1/(right-left)
  else
    return 0
  end
end

function semi_cauchy_variate(u::Float64, left::Float64, right::Float64)
  # two cases, bounded on the left or on the right
  if left != -Inf
    return abs(s_cauchy_variate(u)) + left
  else
    return -abs(s_cauchy_variate(u)) + right
  end
end

function semi_cauchy_pdf(x::Float64, left::Float64, right::Float64)
  if left != -Inf
    if (left <= x)
      return 2*s_cauchy_pdf(x - left)
    else
      return 0
    end
  else
    if (x <= right)
      return 2*s_cauchy_pdf(x - right)
    else
      return 0
    end
  end
end

function proportion(domain, samples)
  filter(x -> (domain[1] < x) && (x < domain[2]), samples) |> length |> x -> x/length(samples)
end


export Sprocket
end
