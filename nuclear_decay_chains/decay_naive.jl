#!/usr/bin/env julia

using DataFrames, Statistics
using Plots, PlotThemes, UnicodePlots
include("decay_helper.jl")


# naive nested-loop O(n^2) atom-counting model
function naive( ; m1=1, m2=1, m3=1, m4=1, disp=false)
  START_ATOMS = 100000      # Bi213
  
  # dependent variables: atom counts
  nBi213 = START_ATOMS
  nPo213 = 0
  nTl209 = 0
  nPb209 = 0
  nBi209 = 0
  # count alpha decay particles
  nAlpha = 0

  h = 1.0       # 1sec time-step width
  tmax = 10000  # total time steps

  # parameters: isotope half-lives (seconds)
  τ1 = 3.253*60*60                # Pb209
  τ1 *= m1
  τ2 = 3.65e-6                    # Po213
  τ2 *= m2
  τ3 = 2.161*60                   # Tl209
  τ3 *= m3
  τ4 = 45.59*60                   # Bi213
  τ4 *= m4

  # decay constants
  λ1 = log(2)/τ1                  # Pb209
  λ2 = log(2)/τ2                  # Po213
  λ3 = log(2)/τ3                  # Tl209
  λ4 = log(2)/τ4                  # Bi213

  # non-branch decay probabilities
  dp1 = 1 - ℯ^(-λ1*h)             # Pb209 -> Bi209
  dp2 = 1 - ℯ^(-λ2*h)             # Po213 -> Pb209
  dp3 = 1 - ℯ^(-λ3*h)             # Tl209 -> Pb209
  dp4 = 1 - ℯ^(-λ4*h)             # Bi213 -> Po213

  # branch decay probabilities
  bp1 = 0.9791                    # Bi213 -> Po213
  bp2 = 0.0209                    # Bi213 -> Tl209

  # atom count vectors
  t_steps = range(1, tmax)
  Bi213 = Vector{Int}()
  Po213 = Vector{Int}()
  Tl209 = Vector{Int}()
  Pb209 = Vector{Int}()
  Bi209 = Vector{Int}()

  # main loop
  for t in t_steps
    push!(Bi213, nBi213)
    push!(Po213, nPo213)
    push!(Tl209, nTl209)
    push!(Pb209, nPb209)
    push!(Bi209, nBi209)

    # compute decay Pb209 -> Bi209
    decay = 0
    for i in range(1, nPb209)
      if rand() <= dp1
        decay += 1
      end
    end
    nPb209 -= decay
    nBi209 += decay

    # compute decay Po213 -> Pb209
    decay = 0
    for i in range(1, nPo213)
      if rand() <= dp2
        decay += 1
      end
    end
    nPo213 -= decay
    nAlpha += decay         # α-particles
    nPb209 += decay

    # compute decay Tl209 -> Pb209
    decay = 0
    for i in range(1, nTl209)
      if rand() <= dp3
        decay += 1
      end
    end
    nTl209 -= decay
    nPb209 += decay

    # compute 2 decay paths
    # Bi213 -> Tl209 &
    # Bi213 -> Po213
    decay = 0
    for i in range(1, nBi213)
      if rand() < dp4
        decay += 1
        if rand() <= bp2    # Bi213 -> Tl209 (rare!)
          nAlpha += 1       # α-particle
          nTl209 += 1
        else                # Bi213 -> Po213 (97.91%)
          nPo213 += 1
        end
      end
    end
    nBi213 -= decay
  end

  # display
  if disp == true
    gr()    # Plots.jl backend
    theme(:ggplot2)

    t = collect(range(start=1, stop=tmax, step=1))
    plt1 = plot(t, Bi213, label="Bi213", xaxis="t", yaxis="atoms", color=:red, linewidth=2)
    plt2 = plot(t, [Po213, Tl209], label=["Po213" "Tl209"], xaxis="t", yaxis="atoms",
                color=[:orange :blue], linewidth=2)
    plt3 = plot(t, [Pb209, Bi209], label=["Pb209" "Bi209"], xaxis="t", yaxis="atoms",
                color=[:green :magenta], linewidth=2)

    #display(plt1)
    #display(plt2)
    #display(plt3)

    # return the isotope vectors
    return [ Bi213, Po213, Tl209, Pb209, Bi209 ]
  end

  # final counts
  final_counts = [ Bi213[end], Po213[end], Tl209[end], Pb209[end], Bi209[end] ]

  # double check sums add up we are using Ints
  if sum(final_counts) != START_ATOMS
    println("uh oh! not ok! vector sum is no good.")
  end

  return final_counts, nAlpha
end # naive()

# batch execute multiple naive() runs 
function batchexec(n = 5; m1=1, m2=1, m3=1, m4=1, disp=false)
  EXECS = n           # total executions

  # results dataframe
  df = DataFrame(Bi213 = Int[], Po213 = Int[], Tl209 = Int[], Pb209 = Int[], Bi209 = Int[])
  # alpha particle vector
  alphav = Vector{Int}() 

  # append naive model results
  for i in range(1, EXECS)
    result, alpha = naive( ; m1, m2, m3, m4)
    push!(df, result)
    push!(alphav, alpha)
  end

  if disp == true
    display(df)
    println()
    #return df
  end

  # summary stats  
  data = DataFrame(Bi213 = Int[], Po213 = Int[], Tl209 = Int[], Pb209 = Int[], Bi209 = Int[])
  # column means
  col_means = [ round(mean(col)) for col in eachcol(df) ]
  push!(data, col_means)
  # column variances
  col_vars = [ round(var(col)) for col in eachcol(df) ]
  push!(data, col_vars)
  # column sd's 
  col_sds  = [ round(std(col)) for col in eachcol(df) ]
  push!(data, col_sds)
  # column max
  col_maxs = [ round(maximum(col)) for col in eachcol(df) ]
  push!(data, col_maxs)
  # column min
  col_mins = [ round(minimum(col)) for col in eachcol(df) ]
  push!(data, col_mins)
  # column spreads
  col_sprds = col_maxs - col_mins
  push!(data, col_sprds)

  return data, round(mean(alphav))
end

# for parameter τ1, Pb209 half-life
function sensitivity_plots_p1()
  adj0 = naive(m1 = 0.95; disp = true)
  unmod = naive(disp = true)
  adj1 = naive(m1 = 1.02; disp = true)
  adj2 = naive(m1 = 1.10; disp = true)
  adj3 = naive(m1 = 1.25; disp = true)

  tmax = 10000
  t = collect(range(start=1, stop=tmax, step=1))

  plt1 = plot(t, [adj0[1], unmod[1], adj1[1], adj2[1], adj3[1]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Bi213 counts, Pb209 half-life (τ1) sensitivity", linewidth=2)
  plt2 = plot(t, [adj0[2], unmod[2], adj1[2], adj2[2], adj3[2]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Po213 counts, Pb209 half-life (τ1) sensitivity", linewidth=2)
  plt3 = plot(t, [adj0[3], unmod[3], adj1[3], adj2[3], adj3[3]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Tl209 counts, Pb209 half-life (τ1) sensitivity", linewidth=2)
  plt4 = plot(t, [adj0[4], unmod[4], adj1[4], adj2[4], adj3[4]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Pb209 counts, Pb209 half-life (τ1) sensitivity", linewidth=2)
  plt5 = plot(t, [adj0[5], unmod[5], adj1[5], adj2[5], adj3[5]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Bi209 counts, Pb209 half-life (τ1) sensitivity", linewidth=2)

  display(plt1)
  display(plt2)
  display(plt3)
  display(plt4)
  display(plt5)
end

# for parameter τ4, Bi213 half-life
function sensitivity_plots_p4()
  adj0 = naive(m4 = 0.95; disp = true)
  unmod = naive(disp = true)
  adj1 = naive(m4 = 1.02; disp = true)
  adj2 = naive(m4 = 1.10; disp = true)
  adj3 = naive(m4 = 1.25; disp = true)

  tmax = 10000
  t = collect(range(start=1, stop=tmax, step=1))
  
  plt1 = plot(t, [adj0[1], unmod[1], adj1[1], adj2[1], adj3[1]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Bi213 counts, Bi213 half-life (τ4)  sensitivity", linewidth=2)
  plt2 = plot(t, [adj0[2], unmod[2], adj1[2], adj2[2], adj3[2]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Po213 counts, Bi213 half-life (τ4) sensitivity", linewidth=2)
  plt3 = plot(t, [adj0[3], unmod[3], adj1[3], adj2[3], adj3[3]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Tl209 counts, Bi213 half-life (τ4) sensitivity", linewidth=2)
  plt4 = plot(t, [adj0[4], unmod[4], adj1[4], adj2[4], adj3[4]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Pb209 counts, Bi213 half-life (τ4) sensitivity", linewidth=2)
  plt5 = plot(t, [adj0[5], unmod[5], adj1[5], adj2[5], adj3[5]],
              label=["-5%" "orig" "+2%" "+10%" "+25%"],legend=:right, xaxis="t", yaxis="atoms",
              title="Bi209 counts, Bi213 half-life (τ4) sensitivity", linewidth=2)
  
  display(plt1)
  display(plt2)
  display(plt3)
  display(plt4)
  display(plt5)
end

# main
clear(5)

# vector similarity
#println(diff_mag(naive, "Pb209", -5))
#println(diff_mag(naive, "Bi213", -5))

#println(diff_mag(naive, "Pb209", 2))
#println(diff_mag(naive, "Bi213", 2))

#println(diff_mag(naive, "Pb209", 10))
#println(diff_mag(naive, "Bi213", 10))

#println(diff_mag(naive, "Pb209", 25))
#println(diff_mag(naive, "Bi214", 25))


# batchexec with original parameters
#orig = batchexec(5, disp=true)
# batchexec but increase param.4 (Bi213 half-life) by 10%
#test1 = batchexec(5, p4 = 1.10, disp=true)
# batchexec but increase param.1 (Pb209 half-life) by 5%
#test2 = batchexec(5, p1 = 1.05, disp=true)

#@show orig
#@show test1
#@show test2
#@show norm(orig-test1)
#@show norm(orig-test2)

#sensitivity_plots_p1()
#sensitivity_plots_p4()
#@time println(batchexec(5, disp=true))
@time println(batchexec(5000, disp=true))