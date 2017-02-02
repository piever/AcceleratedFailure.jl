function test_integration()
    ϕc = [log(10.)]
    t0 = 7.
    t1 = 18.
    pdist = Survival.PGamma(ϕc)
    coefs = Survival.clenshaw_coefs(pdist,t->Survival.ds_dϕ(pdist,t)[1])
    return t0,t1,coefs,pdist
end

t0,t1,coefs,pdist = test_integration()
tempo1 = @benchmark test_integration()

c0 = cdf(pdist,t0)
c1 = cdf(pdist,t1)
r1 = (Survival.clenshaw_asin(c1,coefs)-Survival.clenshaw_asin(c0,coefs))/(c1-c0)
r2 = quadgk(t -> pdf(pdist,t)*Survival.ds_dϕ(pdist,t)[1], t0,t1)[1]/(c1-c0)

tempo2 = @benchmark Survival.clenshaw_asin(rand(), $coefs)
