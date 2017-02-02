using Survival
using Calculus
using Distributions
using Base.Test
using DataFrames
using CSV
using JLD
using BenchmarkTools
using FixedSizeArrays
using Optim

# Check Gamma log likelihood
include("basic_fit.jl")
@test_approx_eq_eps exp(res.minimizer) 10 1e-1
println("Basic Gamma Fit works")

# Check efficient integration, both accuracy and speed
include("integrator.jl")
@test_approx_eq_eps r1 r2 1e-6
println("Efficient integration is fine: $r1 ~ $r2 ")
println("Time elapsed in computing function for integral:")
println(median(tempo1.times)*1e-9)
println("Time elapsed in evaluating function for integral:")
println(median(tempo2.times)*1e-9)

#test cox
include("cox_rossi.jl")
@test_approx_eq_eps hcat(outcome.coefmat.cols[1:3]...) vcat(expected_coefmat.cols[1:3]...) 1e-6
println("Cox regression is fine")
