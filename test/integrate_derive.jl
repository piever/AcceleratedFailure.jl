function test_integration()
    ϕc = [log(10.)]
    t0 = 0.7
    t1 = 1.2
    pdist = Survival.PGamma(ϕc)
    coefs = Survival.clenshaw_coefs(pdist,log)
    return t0,t1,coefs,pdist
end

t0,t1,coefs,pdist = test_integration()
tempo1 = @benchmark test_integration()
println("Time elapsed in computing function for integral:")
println(median(tempo1.times)*1e-9)

c0 = cdf(pdist,t0)
c1 = cdf(pdist,t1)
r1 = (Survival.clenshaw_asin(c1,coefs)-Survival.clenshaw_asin(c0,coefs))/(c1-c0)
r2 = quadgk(t -> pdf(pdist,t)*log(t), t0,t1)[1]/(c1-c0)
@test_approx_eq_eps r1 r2 1e-6
println("Efficient integration is fine: $r1 ~ $r2 ")


tempo2 = @benchmark Survival.clenshaw_asin(rand(), $coefs)
println("Time elapsed in evaluating function for integral:")
println(median(tempo2.times)*1e-9)

ders = Survival.Derivatives(1)
int_coefs = Survival.IntCoefs(pdist)
for s in [Survival.EventWindow(1.),Survival.EventWindow(0.7,Inf),Survival.EventWindow(0.7,1.2)]
    c = 1.
    gradnum = Calculus.gradient(v->Survival.compute_loglik!(ders, s,Survival.PGamma([v[1]]), v[2]),[pdist.params[1],c])
    hesnum = Calculus.hessian(v->Survival.compute_loglik!(ders, s, Survival.PGamma([v[1]]), v[2]),[pdist.params[1],c])
    Survival.compute_loglik!(ders, s,pdist, c)
    Survival.compute_ders!(ders, s, pdist, c, int_coefs )
    @test_approx_eq_eps ders.gradlog gradnum 1e-3
    @test_approx_eq_eps ders.heslog hesnum 1e-3
end
println("Derivation works fine!")
println("Derivatives are $ders")
