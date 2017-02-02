function test_integration()
    ϕc = [log(10.)]
    t0 = 7.
    t1 = 12.
    pdist = Survival.PGamma(ϕc)
    coefs = Survival.clenshaw_coefs(pdist,t->Survival.ds_dϕ(pdist,t,1))
    return t0,t1,coefs,pdist
end

t0,t1,coefs,pdist = test_integration()
tempo1 = @benchmark test_integration()
println("Time elapsed in computing function for integral:")
println(median(tempo1.times)*1e-9)

c0 = cdf(pdist,t0)
c1 = cdf(pdist,t1)
r1 = (Survival.clenshaw_asin(c1,coefs)-Survival.clenshaw_asin(c0,coefs))/(c1-c0)
r2 = quadgk(t -> pdf(pdist,t)*Survival.ds_dϕ(pdist,t, 1), t0,t1)[1]/(c1-c0)
@test_approx_eq_eps r1 r2 1e-6
println("Efficient integration is fine: $r1 ~ $r2 ")


tempo2 = @benchmark Survival.clenshaw_asin(rand(), $coefs)
println("Time elapsed in evaluating function for integral:")
println(median(tempo2.times)*1e-9)

ders = Survival.Derivatives(1)
int_coefs = Survival.IntCoefs(pdist)
for s in [Survival.SurvWindow(10.),Survival.SurvWindow(7.,Inf),Survival.SurvWindow(7.,12.)]
    c = 1.
    Survival.compute_ders!(ders, s, pdist, c, int_coefs )
    value = Survival.compute_loglik(s,pdist, c)
    grad = Calculus.gradient(v->Survival.compute_loglik(s,Survival.PGamma([v[1]]), v[2]),[pdist.params[1],c])
    hes = Calculus.hessian(v->Survival.compute_loglik(s,Survival.PGamma([v[1]]), v[2]),[pdist.params[1],c])
    @test_approx_eq_eps ders.loglik value 1e-3
    @test_approx_eq_eps ders.ds_dϕ[1] grad[1] 1e-3
    @test_approx_eq_eps ders.ds_dc grad[2] 1e-3
    @test_approx_eq_eps ders.d²s_dϕ²[1,1] hes[1,1] 1e-3
    @test_approx_eq_eps ders.d²s_dcdϕ[1] hes[1,2] 1e-3
    @test_approx_eq_eps ders.d²s_dc² hes[2,2] 1e-3
end
println("Derivation works fine!")
println("Derivatives are $ders")
