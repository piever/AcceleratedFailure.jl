# Survival

[![Build Status](https://travis-ci.org/piever/Survival.jl.svg?branch=master)](https://travis-ci.org/piever/Survival.jl)

[![Coverage Status](https://coveralls.io/repos/piever/Survival.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/piever/Survival.jl?branch=master)

[![codecov.io](http://codecov.io/github/piever/Survival.jl/coverage.svg?branch=master)](http://codecov.io/github/piever/Survival.jl?branch=master)

## Basic usage:
Load relevant packages:

```julia
import Survival
using DataFrames
using CSV
```

Load dataset and create event column. Event.time is time, whereas Event.censored is true if the data is censored and false otherwise.

```julia
filepath = joinpath(Pkg.dir("Survival", "test"), "rossi.csv")
rossi = CSV.read(filepath)
rossi[:event] = Survival.Event.(convert(Array, rossi[:week]),convert(Array,rossi[:arrest]) .== 0)
```

Run Cox regression
```julia
outcome = Survival.coxph(event ~ fin+age+race+wexp+mar+paro+prio,rossi)
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
Survival.kaplan_meier(convert(Array,rossi[:event]))
```

or Nelson Aalen estimator for the cumulative hazard:

```julia
Survival.nelson_aalen(convert(Array,rossi[:event]))
```

You can also get the baseline cumulative hazard from the outcome of a Cox regression:

```julia
Survival.nelson_aalen(outcome)
```

## To do list:
- Speeding up Cox with preallocation
- Giving a more extensive output in Cox regression
- Documentation + tests
- Accelerated failure time models!
