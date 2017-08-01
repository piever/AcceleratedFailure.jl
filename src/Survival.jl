module Survival
using DataFrames
using GLM
using PositiveFactorizations
using Distributions
import Base.show
using StatsBase: StatisticalModel, RegressionModel
using ApproxFun
using StaticArrays
import Distributions.pdf, Distributions.logpdf, Distributions.cdf,
       Distributions.quantile, Distributions.rand
import GLM.coef, GLM.predict

include("typedefs.jl")
include("distributions.jl")
include("cox_utils.jl")
include("kaplanmeier.jl")
include("nelsonaalen.jl")
include("optimization.jl")
include("coxph.jl")
include("aft_utils.jl")
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
export coef_reg
export coef_dist

end
