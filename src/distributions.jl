const gamma_pdist = Pdistribution(
(ϕ,t)->pdf(Gamma(ϕ,1),t),
(ϕ,t)->cdf(Gamma(ϕ,1),t),
1,
(ϕ,t)-> [log(t)-polygamma(0,ϕ), (ϕ-1)/t-1],
(ϕ,t)-> [-polygamma(1,ϕ) 1/t; 1/t (ϕ-1)/t^2],
(ϕ,t) -> (quantile(Gamma(ϕ,1.),t)))

function build_auxiliary(pdist::Pdistribution, f::Function, ϕ; degree = 50)
    f1 = Fun(x-> 1/2*cos(x)*f(pdist.quantile(ϕ, (sin(x)+1)/2)), -π/2..π/2, degree)
    g1 = cumsum(f1)
    return g1
end

function dist_integrate(pdist::Pdistribution, f::Function, ϕ, t0, t1; degree = 50)
    g1 = build_auxiliary(pdist, f, ϕ; degree = degree)
    return dist_integrate(pdist, f, g1, ϕ, t0, t1)
end

function dist_integrate(pdist::Pdistribution, f::Function, g::ApproxFun.Fun, ϕ, t0, t1)
    c1 = pdist.cdf(ϕ,t1)
    c0 = pdist.cdf(ϕ,t0)
    return (extrapolate(g,asin(2*c1-1) )-extrapolate(g,asin(2*c0-1)))/(c1-c0)
end

function dist_integrate(pdist::Pdistribution, f::Function, g::ApproxFun.Fun, ϕ, t0)
    c0 = pdist.cdf(ϕ,t0)
    return -extrapolate(g,asin(2*c0-1))/(1-c0)
end
