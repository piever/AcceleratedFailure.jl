using Survival
using DataFrames
using CSV
using BenchmarkTools
using Plots
Plots.plotlyjs(bottom_margin = 5mm, grid = false)
plot(rand(100))
gui()

filepath = joinpath(Pkg.dir("Survival", "examples"), "rossi.csv")
rossi = CSV.read(filepath);
rossi[:event] = Event.(convert(Array, rossi[:week]),convert(Array,rossi[:arrest]) .== 0)

outcome = coxph(event ~ fin+age+race+wexp+mar+paro+prio,rossi)

kaplan_meier(convert(Array,rossi[:event]))

s = Plots.plot(kaplan_meier(convert(Array,rossi[:event]))..., line = :step, label = "kaplan meier")
x,y = nelson_aalen(convert(Array,rossi[:event]))
Plots.plot!(chaz2cdf(x,y)..., line = :step, label = "via nelson aalen")
gui()

x,chaz = nelson_aalen(outcome)

bw = 0.5
Plots.plot(chaz2haz(x,chaz,bw)...)
Plots.plot!(chaz2haz(x,chaz)...)
gui()
