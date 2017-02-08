module Survival
using DataFrames
using GLM
using PositiveFactorizations
using Distributions
import Base.show
using StatsBase: StatisticalModel, RegressionModel
using ApproxFun
using FixedSizeArrays
import Distributions.pdf, Distributions.cdf, Distributions.quantile

include("typedefs.jl")
include("distributions.jl")
include("utils.jl")
include("kaplanmeier.jl")
include("nelsonaalen.jl")
include("optimization.jl")
include("coxph.jl")
include("aft.jl")
include("fast_integral.jl")


export Event
export EventWindow
export coxph
export nelson_aalen
export kaplan_meier
export chaz2cdf
export chaz2haz
export aft
export PGamma

end
