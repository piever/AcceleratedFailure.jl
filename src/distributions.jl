const gamma_pdist = Pdistribution(
(ϕ,t)->pdf(Gamma(ϕ,1),t),
(ϕ,t)->cdf(Gamma(ϕ,1),t),
1,
(ϕ,t)-> [log(t)-polygamma(0,ϕ), (ϕ-1)/t-1],
(ϕ,t)-> [-polygamma(1,ϕ) 1/t; 1/t (ϕ-1)/t^2],
(ϕ,t) -> (quantile(Gamma(ϕ,1.),t)))

function clenshaw_coefs(pdist::Pdistribution, f::Function, ϕ; degree = 50)
    f1 = Fun(x-> 1/2*cos(x)*f(pdist.quantile(ϕ, (sin(x)+1)/2)), -π/2..π/2, degree)
    g1 = cumsum(f1)
    u =  StaticArrays.SVector{length(g1.coefficients)}(g1.coefficients)
    return u
end

clenshaw_corrected(x,c) = clenshaw_halved(4*(asin(2x-1))/π, c)
@generated function clenshaw_halved{N, R<:Real}(x::R, c::StaticArrays.SVector{N,R})
    a,b = :(zero(R)),:(zero(R))
    as = []
    for k = N:-1:2
        ak = Symbol("a",k)
        push!(as, :($ak = $a))
        a = :(muladd(x,$a,c[$k]-$b))
        b = :($ak)
    end
    Expr(:block,
    as...,
    :(muladd(x/2,$a,c[1]-$b)))
end

# function dist_integrate(pdist::Pdistribution, f::Function, ϕ, t0, t1; degree = 50)
#     g2 = build_auxiliary(pdist, f, ϕ; degree = degree)
#     c0 = pdist.cdf(ϕ,t0)
#     c1 = pdist.cdf(ϕ,t1)
#     return (g2(c1) -g2(c0))/(c1-c0)
# end


# function cdfgradhessianvec(pdist,ϕ,ts)
#     f1 = [Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*pdist.gradlog(ϕ,aux(x))[i], aux_inverse(0)..aux_inverse(Inf))
#     for i in 1:pdist.M]
#     g1 = cumsum.(f1)
#     f2 = [Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*
#     (pdist.gradlog(ϕ,aux(x))[i]*pdist.gradlog(ϕ,aux(x))[j]+pdist.heslog(ϕ,aux(x))[i,j]),
#     aux_inverse(0)..aux_inverse(Inf)) for i in 1:pdist.M, j in 1:pdist.M]
#     g2 = cumsum.(f2)
#     grad = zeros(pdist.M+1)
#     hes = zeros(pdist.M+1,pdist.M+1)
#     for t in ts
#         for i in 1:pdist.M
#             grad[i] = (g1[i](aux_inverse(Inf))-g1[i](aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#         end
#         grad[pdist.M+1] = -pdist.pdf(ϕ,t)/(1-pdist.cdf(ϕ,t))
#         hes = -grad*grad'
#         for i in 1:pdist.M
#             for j in 1:pdist.M
#                 hes[i,j] += (g2[i,j](aux_inverse(Inf))-g2[i,j](aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#             end
#         end
#         hes[end, :] += grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#         hes[:,end] = hes[end, :]
#     end
#     return grad, hes
# end
