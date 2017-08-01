# Indefinite integral of pdist.pdf*f is clenshaw_asin(cdf(x),clenshaw_coefs)

function clenshaw_coefs{N}(pdist::Distribution, f::Function, degreetype::Val{N} = Val{50}())
    f1 = Fun(x-> 1/2*cos(x)*f(quantile(pdist, (sin(x)+1)/2)), -π/2..π/2, N-1)
    g1 = cumsum(f1)
    len = min(length(g1.coefficients),N)
    (len < length(g1.coefficients)) && warn("I'm truncating the coefficients!")
    u =  SVector{N, Float64}(vcat(g1.coefficients[1:len],zeros(N-len)))
    return u
end

#clenshaw_asin(x,c) = ApproxFun.clenshaw(2*(asin(2x-1))/π, c)
clenshaw_asin(x,c) = clenshaw_halved(4*(asin(2x-1))/π, c)
@generated function clenshaw_halved{N, R<:Real}(x::R, c::SVector{N})
    Expr(:block,
    :((a,b) = (zero(R), zero(R))),
    (:((a,b) = (muladd(x,a,c[$k]-b),a)) for k = N:-1:2)...,
    :(muladd(x/2,a,c[1]-b)))
end
