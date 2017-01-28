ϕc = 10.
t0 = 7.
t1 = 18.

r1 = Survival.dist_integrate(Survival.gamma_pdist,t -> Survival.gamma_pdist.gradlog(ϕc,t)[1], ϕc, t0, t1)
r2 = quadgk(t -> Survival.gamma_pdist.pdf(ϕc,t)*Survival.gamma_pdist.gradlog(ϕc,t)[1], t0,t1)[1]/
(Survival.gamma_pdist.cdf(ϕc,t1)-Survival.gamma_pdist.cdf(ϕc,t0))
