filepath = joinpath(Pkg.dir("Survival", "examples"), "rossi.csv")
rossi = CSV.read(filepath; nullable = false)
rossi[:event] = Event.(rossi[:week],rossi[:arrest] .== 0)

outcome = coxph(event ~ fin+age+race+wexp+mar+paro+prio,rossi; tol = 1e-8)

expected_coefmat = JLD.load("expected_coefmat.jld", "expected_coefmat")
