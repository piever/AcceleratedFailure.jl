using Survival
using Calculus
using Distributions
using Base.Test
using DataFrames
using CSV
using JLD
using BenchmarkTools
using FixedSizeArrays

# Check aft
include("aft_test.jl")

# Check efficient integration, both accuracy and speed
#include("integrate_derive.jl")


#test cox
include("cox_rossi.jl")
@test_approx_eq_eps hcat(outcome.coefmat.cols[1:3]...) vcat(expected_coefmat.cols[1:3]...) 1e-6
println("Cox regression is fine")
