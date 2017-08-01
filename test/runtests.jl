using Survival
using Distributions
using Base.Test
using DataFrames
using JLD

# Check aft
include("aft_test.jl")

#test cox
include("cox_rossi.jl")
