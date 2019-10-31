include("../src/Sprocket.jl")
include("../src/Exact.jl")
using .Sprocket
using .Exact
using JuMP
using SpecialFunctions
using Random


function test_one()
  function build_problem()
    #create an empty problem
    problem = Sprocket.Problem()

    # add one real valued random variable to the problem
    xi = Sprocket.Variable(name=:xi, size=(),type = Sprocket.Random())

    Sprocket.add_variable(problem,xi)

    # add another
    lam = Sprocket.Variable(name=:lam, size=(),type = Sprocket.Random())
    Sprocket.add_variable(problem,lam)

    #set the domain of the random variable
    Sprocket.set_domain(problem,(Sprocket.Point(Dict(xi =>0.0,lam => 0.0 )), Sprocket.Point(Dict(xi => 1.0, lam => 1.0))))

    # objective function is defined through the cutting plane oracle
    function my_cutting_plane_oracle()
      function oracle(vars)
        value = vars[:xi]^2 + 3*vars[:lam]^2 + 0.1*vars[:lam]*vars[:xi]
        grad = deepcopy(vars)
        grad[:xi] = 2*vars[:xi] + 0.1*vars[:lam]
        grad[:lam] = 0.1*vars[:xi] + 6*vars[:lam]
        return (value,grad)
      end

      return oracle
    end

    # distribution of random variables is defined through the measure oracle
    function my_measure_oracle()

      function oracle(vars_1,vars_2)
        left_xi = vars_1[:xi]
        right_xi = vars_2[:xi]
        if vars_1[:xi] <= 0
          left_xi = 0.0
        end
        if vars_1[:xi] >= 1
          left_xi = 1.0
        end
        if vars_2[:xi] <= 0
          right_xi = 0.0
        end
        if vars_2[:xi] >= 1
          right_xi = 1.0
        end
        factor_xi = right_xi - left_xi

        left_lam = vars_1[:lam]
        right_lam = vars_2[:lam]
        if vars_1[:lam] <= 0
          left_lam = 0.0
        end
        if vars_1[:lam] >= 1
          left_lam = 1.0
        end
        if vars_2[:lam] <= 0
          right_lam = 0.0
        end
        if vars_2[:lam] >= 1
          right_lam = 1.0
        end
        factor_lam = right_lam - left_lam

        fresh = deepcopy(vars_1)
        fresh[:xi] = (right_xi^2 - left_xi^2)*factor_lam/2
        fresh[:lam] = (right_lam^2 - left_lam^2)*factor_xi/2
        return (factor_lam*factor_xi,fresh)
      end

    return oracle
    end

      problem.c_oracle = my_cutting_plane_oracle()
      problem.m_oracle = my_measure_oracle()

    return problem
  end

  p = build_problem()
  samples = []

  for i = 1:3000
    push!(samples,Sprocket.random_sample(p.m_oracle,100,p.domain))
  end

  reduce((x,y) -> x + y, samples) / length(samples)

  Sprocket.proportion((0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars)), samples)
  p.m_oracle(0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars))[1]

  println(Sprocket.proportion((0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars)), samples))
  println(p.m_oracle(0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars))[1])
end


function test_two()
  Random.seed!(0)
  function build_problem()
    #create an empty problem
    problem = Sprocket.Problem()

    # add one real valued random variable to the problem
    xi = Sprocket.Variable(name=:xi, size=(),type=Sprocket.Random())
    Sprocket.add_variable(problem,xi)

    #set the domain of the random variable
    # Sprocket.set_domain(problem,Exact.rect_hull,(Sprocket.Point(Dict(xi =>-Inf)), Sprocket.Point(Dict(xi => Inf))))
    Sprocket.set_domain(problem,(Sprocket.Point(Dict(xi =>-Inf)), Sprocket.Point(Dict(xi => Inf))))

    # objective function is defined through the cutting plane oracle
    # This function has unbounded support but is also Lipschitz continuous with Lipschitz constant 1
    function my_cutting_plane_oracle()
      function oracle(vars)
        value = -log(2) - vars[:xi]/2 + log(exp(vars[:xi]) + 1)
        grad = deepcopy(vars)
        grad[:xi] = 1/2 - 1/(1 + exp(vars[:xi]))
        return (value,grad)
      end
    return oracle
    end

    # distribution of random variables is defined through the measure oracle
    # below is the oracle for the standard normal distribution
    function my_measure_oracle()
      function oracle(vars_1,vars_2)
        probability = 0.5*(erf(vars_2[:xi]/sqrt(2)) - erf(vars_1[:xi]/sqrt(2)))
        first_order =  deepcopy(vars_1)
        first_order[:xi] =  (exp(-vars_1[:xi]^2/2) - exp(-vars_2[:xi]^2/2))/sqrt(2*pi)

        return (probability,first_order)
      end

    return oracle
    end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

    return problem
  end

  p = build_problem()
  samples = []

  for i = 1:5000
    push!(samples,Sprocket.random_sample(p.m_oracle,6,p.domain))
  end

  reduce((x,y) -> x + y, samples) / length(samples)

  Sprocket.proportion((-0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars)), samples)
  p.m_oracle(-0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars))[1]

  println(Sprocket.proportion((-0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars)), samples))
  println(p.m_oracle(-0.5*Sprocket.unit(p.vars), 0.7*Sprocket.unit(p.vars))[1])

end


function test_three()
  function build_problem()
    #create an empty problem
    problem = Sprocket.Problem()

    # add one real valued random variable to the problem
    xi = Sprocket.Variable(name=:xi, size=(),type=Sprocket.Random())
    Sprocket.add_variable(problem,xi)

    #set the domain of the random variable
    # Sprocket.set_domain(problem,Exact.rect_hull,(Sprocket.Point(Dict(xi =>-Inf)), Sprocket.Point(Dict(xi => Inf))))
    Sprocket.set_domain(problem,(Sprocket.Point(Dict(xi => 0.0)), Sprocket.Point(Dict(xi => Inf))))

    # objective function is defined through the cutting plane oracle
    # This function has unbounded support but is also Lipschitz continuous with Lipschitz constant 1
    function my_cutting_plane_oracle()
      function oracle(vars)
        value = -log(2) - vars[:xi]/2 + log(exp(vars[:xi]) + 1)
        grad = deepcopy(vars)
        grad[:xi] = 1/2 - 1/(1 + exp(vars[:xi]))
        return (value,grad)
      end
    return oracle
    end

    # distribution of random variables is defined through the measure oracle
    # below is the oracle for the standard normal distribution
    function my_measure_oracle()
      function oracle(vars_1,vars_2)
        left = vars_1[:xi]
        if vars_1[:xi] <= 0
          left = 0.0
        end

        probability = (erf(vars_2[:xi]/sqrt(2)) - erf(left/sqrt(2)))
        first_order =  deepcopy(vars_1)
        first_order[:xi] =  2*(exp(-vars_1[:xi]^2/2) - exp(-left^2/2))/sqrt(2*pi)

        return (probability,first_order)
      end

    return oracle
    end

    problem.c_oracle = my_cutting_plane_oracle()
    problem.m_oracle = my_measure_oracle()

    return problem
  end

  p = build_problem()
  samples = []

  for i = 1:1000
    push!(samples,Sprocket.random_sample(p.m_oracle,10,p.domain))
  end

  reduce((x,y) -> x + y, samples) / length(samples)

  Sprocket.proportion((1.5*Sprocket.unit(p.vars), 2.0*Sprocket.unit(p.vars)), samples)
  p.m_oracle(1.5*Sprocket.unit(p.vars), 2.0*Sprocket.unit(p.vars))[1]

  println( Sprocket.proportion((1.5*Sprocket.unit(p.vars), 2.0*Sprocket.unit(p.vars)), samples) )

  println( p.m_oracle(1.5*Sprocket.unit(p.vars), 2.0*Sprocket.unit(p.vars))[1] )


end
