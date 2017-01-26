const gamma = Pdistribution(
(ϕ,t)->pdf(Gamma(ϕ,1),t),
(ϕ,t)->cdf(Gamma(ϕ,1),t),
1,
(ϕ,t)-> [log(t)-polygamma(0,ϕ), (ϕ-1)/t-1],
(ϕ,t)-> [-polygamma(1,ϕ) 1/t; 1/t (ϕ-1)/t^2])
