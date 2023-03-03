#!/usr/bin/env julia

using Plots, PlotThemes, UnicodePlots
include("decay_helper.jl")


# markov chain model
# all args are optional: m1, m2, etc. are parameter multipliers
function markov( ; m1=1, m2=1, m3=1, m4=1, disp=false)
  START_ATOMS = 100000  # Bi213 initial atoms
  tmax = 10000          # total time steps
  h = 1.0               # 1s time-step width

  # isotope half-lives (seconds)
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
  dp1 = 1 - ℯ^(-λ1*h)             # Pb209
  dp2 = 1 - ℯ^(-λ2*h)             # Po213
  dp3 = 1 - ℯ^(-λ3*h)             # Tl209
  dp4 = 1 - ℯ^(-λ4*h)             # Bi213

  # branch decay probabilities
  bp1 = 0.9791                    # Bi213 -> Po213
  bp2 = 0.0209                    # Bi213 -> Tl209

  # joint probabilities
  # ASSUMPTION: independent events!
  dp5 = dp4 * bp1                 # Bi213 -> Po213
  dp6 = dp4 * bp2                 # Bi213 -> Tl209

  # 5x5 markov chain transition matrix
  #     Bi213,Po213,Tl209,Pb209,Bi209
  A = [ 1-dp4     0     0     0     0;
          dp5 1-dp2     0     0     0;
          dp6     0 1-dp3     0     0;
            0   dp2   dp3 1-dp1     0;
            0     0     0   dp1     1]

  # initial; matrix exponentiation
  x0 = [START_ATOMS; 0; 0; 0; 0;]
  results = A^tmax * x0

  # display
  if disp == true
    gr()
    theme(:ggplot2)
    out = zeros((tmax, length(x0)))
    out[1, :] = x0
    for n in range(start=1, stop=tmax, step=1)
      out[n, :] = A*x0
      x0 = out[n, :]
    end
    t = collect(range(start=1, stop=tmax, step=1))

    # UnicodePlots
    #display(lineplot(t, out[:,1], title="Bismuth-213", color=:red))
    #display(lineplot(t, out[:,2], title="Polonium-213", color=:green))
    #display(lineplot(t, out[:,3], title="Thallium-209", color=:green))
    #display(lineplot(t, out[:,4], title="Lead-209", color=:blue))
    #display(lineplot(t, out[:,5], title="Bismuth-209", color=:blue))
    
    plt1 = plot(t, out[:,1], label="Bi213", color=:red,
           xaxis="t", yaxis="atoms", linewidth=2)
    plt2 = plot(t, [out[:,2], out[:,3]], label=["Po213" "Tl209"], color=[:orange :blue],
           xaxis="t", yaxis="atoms", linewidth=2)
    plt3 = plot(t, [out[:,4], out[:,5]], label=["Pb209" "Bi209"], color=[:green :magenta],
           xaxis="t", yaxis="atoms", linewidth=2)
    
    display(plt1)
    display(plt2)
    display(plt3)

    return out
  end
  
  return [ round(item) for item in results ]
end # markov()

# increase a parameter by a specified percentage
# compute similarity between unmodified and adjusted final step vectors
# vectors represent distribution of atoms in 5 sets: Bi213 Po213 Tl209 Pb209 Bi209
# diff_mag() returns magnitude (norm) of difference between unmod. & adj. vectors
function vector_similarity()
  println(diff_mag(markov, "Pb209", -5))
  #println(diff_mag(markov, "Po213", 2))
  #println(diff_mag(markov, "Tl209", 2))
  println(diff_mag(markov, "Bi213", -5))
 
  println(diff_mag(markov, "Pb209", 2))
  #println(diff_mag(markov, "Po213", 2))
  #println(diff_mag(markov, "Tl209", 2))
  println(diff_mag(markov, "Bi213", 2))

  println(diff_mag(markov, "Pb209", 10))
  #println(diff_mag(markov, "Po213", 10))
  #println(diff_mag(markov, "Tl209", 10))
  println(diff_mag(markov, "Bi213", 10))

  println(diff_mag(markov, "Pb209", 25))
  #println(diff_mag(markov, "Po213", 25))
  #println(diff_mag(markov, "Tl209", 25))
  println(diff_mag(markov, "Bi213", 25))
end


# adjust a parameter by a specified multiplier
# plot multiple data series on the same graph
# param.1 = τ1 = Pb209 half-life
function sensitivity_plots_τ1()
  out0 = markov(m1 = 0.95, disp = true)
  out1 = markov(disp = true)
  out2 = markov(m1 = 1.02, disp = true)
  out3 = markov(m1 = 1.10, disp = true)
  out4 = markov(m1 = 1.25, disp = true)

  t = collect(range(start=1, stop=10000, step=1))
  plt1 = plot(t, [out0[:,1], out1[:,1], out2[:,1], out3[:,1], out4[:,1]],
              title="Bi213 count, Pb209 half-life (τ1) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt2 = plot(t, [out0[:,2], out1[:,2], out2[:,2], out3[:,2], out4[:,2]],
              title="Po213 count, Pb209 half-life (τ1) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt3 = plot(t, [out0[:,3], out1[:,3], out2[:,3], out3[:,3], out4[:,3]],
              title="Tl209 count, Pb209 half-life (τ1) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt4 = plot(t, [out0[:,4], out1[:,4], out2[:,4], out3[:,4], out4[:,4]],
              title="Pb209 count, Pb209 half-life (τ1) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt5 = plot(t, [out0[:,5], out1[:,5], out2[:,5], out3[:,5], out4[:,5]],
              title="Bi209 count, Pb209 half-life (τ1) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  display(plt1)
  display(plt2)
  display(plt3)
  display(plt4)
  display(plt5)
end

# adjust a parameter by a specified multiplier
# plot multiple data series on the same graph
# param.4 = τ4 = Bi213 half-life
function sensitivity_plots_τ4()
  out0 = markov(m4 = 0.95, disp = true)
  out1 = markov(disp = true)
  out2 = markov(m4 = 1.02, disp = true)
  out3 = markov(m4 = 1.10, disp = true)
  out4 = markov(m4 = 1.25, disp = true)

  t = collect(range(start=1, stop=10000, step=1))
  plt1 = plot(t, [out0[:,1], out1[:,1], out2[:,1], out3[:,1], out4[:,1]],
              title="Bi213 count, Bi213 half-life (τ4) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)
  
  plt2 = plot(t, [out0[:,2], out1[:,2], out2[:,2], out3[:,2], out4[:,2]],
              title="Po213 count, Bi213 half-life (τ4) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)
  
  plt3 = plot(t, [out0[:,3], out1[:,3], out2[:,3], out3[:,3], out4[:,3]],
              title="Tl209 count, Bi213 half-life (τ4) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt4 = plot(t, [out0[:,4], out1[:,4], out2[:,4], out3[:,4], out4[:,4]],
              title="Pb209 count, Bi213 half-life (τ4) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  plt5 = plot(t, [out0[:,5], out1[:,5], out2[:,5], out3[:,5], out4[:,5]],
              title="Bi209 count, Bi213 half-life (τ4) sensitivity", xaxis="t", yaxis="atoms",
              label = ["-5%" "orig" "+2%" "+10%" "+25%"], legend=:right)

  display(plt1)
  display(plt2)
  display(plt3)
  display(plt4)
  display(plt5)
end

clear(32)
vector_similarity()
#sensitivity_plots_τ1()
#sensitivity_plots_τ4()
#markov(disp = true)