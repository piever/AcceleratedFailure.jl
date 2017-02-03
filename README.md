# Survival

[![Build Status](https://travis-ci.org/piever/Survival.jl.svg?branch=master)](https://travis-ci.org/piever/Survival.jl)

[![Coverage Status](https://coveralls.io/repos/piever/Survival.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/piever/Survival.jl?branch=master)

[![codecov.io](http://codecov.io/github/piever/Survival.jl/coverage.svg?branch=master)](http://codecov.io/github/piever/Survival.jl?branch=master)

## Setting up
To install the package, simply run
```julia
Pkg.clone("https://github.com/piever/Survival.jl.git")
```
in the Julia REPL.

## Basic usage:
Load relevant packages:

```julia
using Survival
using DataFrames
using CSV
```

Load dataset and create event column. Event.time is time, whereas Event.censored is true if the data is censored and false otherwise.

```julia
filepath = joinpath(Pkg.dir("Survival", "examples"), "rossi.csv")
rossi = CSV.read(filepath; nullable = false)
rossi[:event] = Event.(rossi[:week], rossi[:arrest].== 0)
```

Run Cox regression
```julia
outcome = coxph(event ~ fin+age+race+wexp+mar+paro+prio,rossi)
```
And you should get this outcome (computed with Efron method for ties):
```
Model: Cox; Formula: event ~ fin + age + race + wexp + mar + paro + prio

        Estimate Std.Error   z value Pr(>|z|)
fin    -0.378636  0.191352  -1.97874   0.0478
age   -0.0570578 0.0219674  -2.59738   0.0094
race    0.314349  0.308045   1.02046   0.3075
wexp   -0.152399  0.212026 -0.718774   0.4723
mar    -0.431883  0.381726   -1.1314   0.2579
paro  -0.0850866  0.195743 -0.434686   0.6638
prio   0.0907363 0.0286157   3.17086   0.0015
```

This package also includes Kaplan Meier estimator of the cumulative function, which takes a vector of events as input:

```julia
x,cumulative = kaplan_meier(rossi[:event])
```
In the output `x` will be the array with the times of death, `cumulative` the estimated cdf.

To visualize it (using Plots.jl) you can simply type:

```julia
using Plots
plot(x,cumulative, line = :step)
```

Nelson Aalen estimator for the cumulative hazard:

```julia
x, chaz = nelson_aalen(rossi[:event])
```

To check that everything went well, you may want to verify that `cumulative = 1 - exp.(-chaz)`

```julia
plot(x,cumulative, line = :step)
plot!(x,1-exp.(-chaz), line = :step)
```

You can also get the baseline cumulative hazard from the outcome of a Cox regression:

```julia
x, chaz = nelson_aalen(outcome)
```
In turn, it's always possible to get the cdf from the cumulative hazard:

```julia
chaz2cdf(x,chaz)
```

and the same is true for the hazard:

```julia
chaz2haz(x,chaz)
chaz2haz(x,chaz, bw)

```

where, in the second case, `bw` is a parameter used for smoothing.

## Accelerated Failure Time models

This package also supports [accelerated failure time models](https://en.wikipedia.org/wiki/Accelerated_failure_time_model).

```julia
using Distributions
using DataFrames
using Survival
```

Let's generate a fake dataset:

```julia
N = 50000
x = randn(N)
y = randn(N)
z = randn(N)
t1 = rand.(Gamma.(10,1*exp.(x-0.3y)))
t2 = rand(Gamma(15,1),N)
W = [(t2[i]>t1[i]) ? Event(t1[i], false) : Event(t2[i], true) for i = 1:N]
df = DataFrame(x = x, y = y, z = z, a = W)
```

Let's specify the formula and distribution (only `PGamma` is implemented so far, where `PGamma(params) = Gamma(params[1],1)`):

```julia
res = aft(a ~ 1 + x +y + x*y+ z, df, PGamma(); tol = 1e-3)
```

The outcome should look something like this:

```
Model: Accelerated Failure Time, dist = Survival.PGamma(params=[2.33543]);
Formula: a ~ 1 + x + y + x * y + z

                Estimate  Std.Error   z value Pr(>|z|)
params1          2.33543  0.0237773   98.2207   <1e-99
(Intercept)   -0.0348241  0.0259545  -1.34174   0.1797
x                1.00741 0.00711357   141.617   <1e-99
y              -0.294449 0.00622455  -47.3045   <1e-99
z            -0.00455039 0.00522966 -0.870112   0.3842
x & y          0.0116439 0.00676265   1.72179   0.0851
```

## Using DataFrames
`nelson_aalen`, `kaplan_meier`, Cox regression and accelerated failure time models can be used with Nullable Arrays. The rows missing relevant data will be discarded.
