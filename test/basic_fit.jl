N = 100000
t1 = rand(Gamma(10,1/10),N)
t2 = rand(Gamma(15,1/12),N)

W = [(t2[i]>t1[i]) ? EventWindow(t1[i], t1[i]) : EventWindow(t2[i], Inf) for i = 1:N];
res = optimize(t -> -sum(Survival.compute_loglik.(W, [PGamma([t])], zeros(N))), log(1.),log(20.))
