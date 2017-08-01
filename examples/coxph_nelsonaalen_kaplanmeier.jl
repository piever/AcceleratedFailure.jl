#Recommended use of the example: run step by step in Juno

using Survival
using DataFrames
using BenchmarkTools
using Plots
gr()

# Load DataFrame and create "event" column
filepath = joinpath(Pkg.dir("Survival", "examples"), "rossi.csv")
rossi = readtable(filepath)
rossi[:event] = Event.(rossi[:week],rossi[:arrest] .== 0)

# kaplan_meier and nelson_aalen method to estimate survival and cumulative hazard respectively
plot(kaplan_meier(rossi[:event])..., line = :step)

plot(nelson_aalen(rossi[:event])..., line = :step)

# Run Cox regression
outcome = coxph(@formula(event ~ fin+age+race+wexp+mar+paro+prio),rossi)

# Test efficiency of the implementation:
using BenchmarkTools
@benchmark outcome = coxph(@formula(event ~ fin+age+race+wexp+mar+paro+prio),rossi)

# The outcome of the regression is stored in outcome.coefmat
outcome.coefmat

#Other relevant quantities are stored in outcome as well:
outcome.loglik #log-likelihood of the fitted parameters (for model comparison)
outcome.score # gradient of log-likelihood: should be close to 0 (safety check)
outcome.fischer_info #fischer information matrix

eigvals(outcome.fischer_info) # Safety check: these eigen-values must be positive

# Get baseline cumulative hazard from Cox regression
x,chaz = nelson_aalen(outcome)

# Use cumulative hazard to get hazard (bw is smoothing parameter)
bw = 0.5
plot(chaz2haz(x,chaz,bw)...)
plot!(chaz2haz(x,chaz)...)

# Use cumulative hazard to get cdf
plot(chaz2cdf(x,chaz), line = :step)
