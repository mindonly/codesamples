#!/usr/bin/env julia

using OrdinaryDiffEq

# area of a circle
function A(u)
  r = u
  return π*r^2
end

# du/dt = f(u,p,t)
# -- ODE right-hand side -- consumed by ODEProblem()
# u = water height in conical tank, ft.
# p = <parameters> default: Null
# t = time step, mins.
function f(u, p, t)
  r = 0.1   # radius bottom hole, ft.
  g = 32.1  # Earth grav. constant, ft.

  # Toricelli's Law - relates fluid speed flowing from an orifice
  return -0.6π*r^2 * sqrt(2g) * (u^(1/2) / A(u))
end

# ODE problem setup
function conical_tank(func)
  u0 = 8.0    # initial water height, ft.
  sols = []

  for t in vcat([30, 60, 120, 240, 480, 960], collect(range(1500, stop=1510)))
    prob = ODEProblem(func, u0, (0.0, t))

    try
      sol = solve(prob, Tsit5(), abstol=1e-9, reltol=1e-6)
      push!(sols, sol)
    catch err
      @warn err, "t = $t"
      break
    end
  end

  return sols
end

# main
function report()
  println("\n")

  solutions = conical_tank(f)
  summaries = [ ( sol.t[end], sol.u[length(sol.u)-3:end] ) for sol in solutions ]
  
  return summaries
end

@time report()
