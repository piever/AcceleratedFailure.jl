# Indefinite integral of pdist.pdf*f is clenshaw_asin(cdf(x),clenshaw_coefs)

function clenshaw_coefs{N}(pdist::Pdistribution, f::Function, ϕ, output_type::Val{N} = Val{50}())
    f1 = Fun(x-> 1/2*cos(x)*f(pdist.quantile(ϕ, (sin(x)+1)/2)), -π/2..π/2, N-1)
    g1 = cumsum(f1)
    len = min(length(g1.coefficients),N)
    (len < length(g1.coefficients)) && warn("I'm truncating the coefficients!")
    u =  Vec{N, Float64}(vcat(g1.coefficients[1:len],zeros(N-len)))
    return u
end

#clenshaw_asin(x,c) = ApproxFun.clenshaw_fast(2*(asin(2x-1))/π, c)
clenshaw_asin(x,c) = clenshaw_halved(4*(asin(2x-1))/π, c)
@generated function clenshaw_halved{N, R<:Real}(x::R, c::Vec{N,R})
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
