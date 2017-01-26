using Survival
using Calculus
using Distributions
using Base.Test
using DataFrames
using CSV
using JLD

# write your own tests here

# Test gradient and hessian in aft_gradhes
include("aft_derivative.jl")
@test_approx_eq_eps grad grad1 1e-6
@test_approx_eq_eps hes hes1 1e-6


#test cox
include("cox_rossi.jl")
@test_approx_eq_eps hcat(outcome.coefmat.cols[1:3]...) vcat(expected_coefmat.cols[1:3]...) 1e-6
