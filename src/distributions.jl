const gamma_pdist = Pdistribution(
(ϕ,t)->pdf(Gamma(ϕ,1),t),
(ϕ,t)->cdf(Gamma(ϕ,1),t),
1,
(ϕ,t)-> [log(t)-polygamma(0,ϕ), (ϕ-1)/t-1],
(ϕ,t)-> [-polygamma(1,ϕ) 1/t; 1/t (ϕ-1)/t^2],
(ϕ,t) -> (quantile(Gamma(ϕ,1.),t)))

# MANCANO::
# get_int_coefs!
# compute_ders!
# compute_loglik
