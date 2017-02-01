N = 100000
t1 = rand(Gamma(10,1),N)
t2 = rand(Gamma(15,1),N)

W = [(t2[i]>t1[i]) ? Survival.SurvWindow(t1[i], true) : Survival.SurvWindow(t2[i], false) for i = 1:N];
res = optimize(t -> -sum(Survival.compute_loglik.(W, [Survival.PGamma([t])], zeros(N))), log(1.),log(20.))
