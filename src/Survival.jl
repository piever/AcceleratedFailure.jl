module Survival
using DataFrames
using GLM
using Optim
using TensorOperations
using PositiveFactorizations
using Distributions
import ForwardDiff
import Base.show
using StatsBase: StatisticalModel, RegressionModel

type EventHistoryModel <: RegressionModel
    model::AbstractString
    formula::Formula
    eventtype::DataType
    #coef::Vector{Float64}
    coefmat::CoefTable
    #IC::Matrix{Float64}
    #invhess::Matrix{Float64}
    #vcov::Matrix{Float64}
    #opt::Vector
    #grad::Vector{Float64}
    #X::Matrix{Float64}
    eventtime::Matrix
    #chaz::Matrix{Float64}
end

# coefnames !!!

function show(io::IO, obj::EventHistoryModel)
    print(io,"\nModel: ", obj.model, "; ", obj.formula,"\n")
    n = size(obj.eventtime,1)
    events::Int = sum(obj.eventtime[:,2])
    print(io,"\nn=",n,", events=",events,"\n\n")
    print(io,obj.coefmat)
end

include("utils.jl")
include("kaplanmeier.jl")
include("parametric.jl")
include("optimization.jl")
include("coxph.jl")

export Event
export coxph

end
