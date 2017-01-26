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

## Using DataFrames
`nelson_aalen`, `kaplan_meier` and Cox regression can be used with Nullable Arrays. The rows missing relevant data will be discarded.

## To do list:
- Documentation + tests + examples
- Accelerated failure time models
