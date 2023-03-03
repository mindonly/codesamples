#!/usr/bin/env julia

using CodeTracking


# clear screen
clear(N) = for i = 1:N
  print("\n")
end

# func = function obj. naive() or markov()
# iso = isotope "Po213" "Pb209" etc.
# percent = integer param. adjustment 
function diff_mag(func, iso, percent)
  if iso == "Pb209"
    orig = func(m1 = 1)
    test = func(m1 = 1 + percent/100)
  elseif iso == "Po213"
    orig = func(m2 = 1)
    test = func(m2 = 1 + percent/100)
  elseif iso == "Tl209"
    orig = func(m3 = 1)
    test = func(m3 = 1 + percent/100)
  elseif iso == "Bi213"
    orig = func(m4 = 1)
    test = func(m4 = 1 + percent/100)
  end

  println("\n$percent% increase, $iso half-life")
  @show orig
  @show test
  println("       [Bi213  Po213  Tl209  Pb209  Bi209]")

  # use [1] for atom counts vector
  # [2] is alpha particle count
  vecdiff = norm(orig[1]-test[1])
  #@show vecdiff

  return vecdiff
end

function element_report(τ)
  τ = τ           # half-life (s)
  λ = log(2)/τ    # decay constant
  h = 1.0         # 1sec step-width
  
  println("probability functions:")
  a(h) = 1 - 2^(-h/τ)   # Newman
  b(h) = 1 - ℯ^(-λ*h)   # Shultis/Faw
  println(@code_string a(0))
  println(@code_string b(0))
  println()
  
  @show τ
  @show λ
  println()
  @show λ - a(1)
  @show λ - b(1)
  println()

  @show a(1)
  @show b(1)
  @show a(100)
  @show b(100)
end

function solcheck(τ)
  τ = τ           # half-life (s)
  λ = log(2)/τ    # decay constant
  N0 = 100000     # initial # of atoms
  tmax = 10000    # total time steps

  println("\ndiffeq solution:")
  N(t) = N0 * ℯ^(-λ*t)
  println(@code_string N(0))
  println("N($tmax) = ", N(tmax))
end

# main
#println("\nPb209\n-----")
#element_report(3.3*60*60)

#println("\nTl209\n-----")
#element_report(2.2*60)

#println("\nBi213\n-----")
#element_report(45.59*60)
#solcheck(45.59*60)

#println("\nRn220\n-----")
#element_report(55.61)
#solcheck(55.61)

#element_report(3.65e-6)