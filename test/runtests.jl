using Survival
using Calculus
using Distributions
using Base.Test
using DataFrames
using CSV
using JLD
using BenchmarkTools

# Check aft
include("aft_test.jl")

#test cox
include("cox_rossi.jl")
