function test_integration()
    ϕc = 10.
    t0 = 7.
    t1 = 18.
    pdist = Survival.gamma_pdist
    coefs = Survival.clenshaw_coefs(Survival.gamma_pdist,t -> Survival.gamma_pdist.gradlog(ϕc,t)[1], ϕc)
    return ϕc,t0,t1,coefs,pdist
end

ϕc,t0,t1,coefs,pdist = test_integration()
tempo1 = @benchmark test_integration()

c0 = pdist.cdf(ϕc,t0)
c1 = pdist.cdf(ϕc,t1)
r1 = (Survival.clenshaw_corrected(c1,coefs)-Survival.clenshaw_corrected(c0,coefs))/(c1-c0)
r2 = quadgk(t -> Survival.gamma_pdist.pdf(ϕc,t)*Survival.gamma_pdist.gradlog(ϕc,t)[1], t0,t1)[1]/
(Survival.gamma_pdist.cdf(ϕc,t1)-Survival.gamma_pdist.cdf(ϕc,t0))


function f(k)
    for i = 1:1000
        Survival.clenshaw_corrected(rand(), k)
    end
end
tempo2 = @benchmark f(coefs)
