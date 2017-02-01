function compute_loglik(s::SurvWindow, dist::Distribution, coef)
    if s.instant
        return log(pdf(dist,s.t₀*exp(-coef)))-coef
    elseif s.t₁ < Inf
        return log(cdf(dist,s.t₁*exp(-coef))-cdf(dist,s.t₀*exp(-coef)))
    else
        return log(1-cdf(dist,s.t₀*exp(-coef)))
    end
end

function compute_ders!(ders, s, pdist::Distribution, c, int_coefs)
    ϕ = pdist.params[1]
    ders.loglik = compute_loglik(s, pdist, c)
    τ₁,τ₀ = [s.t₁,s.t₀]*exp(-c)
    if s.instant
        ders.ds_dϕ[:]    = ds_dϕ(pdist, τ₀)
        ders.ds_dc       = ds_dc(pdist, τ₀)
        ders.ds²dϕ²[:,:] = ds²_dϕ²(pdist, τ₀)
        ders.ds²_dcdϕ[:] = ds²_dcdϕ(pdist, τ₀)
        ders.ds²_dc²     = ds²_dc²(pdist, τ₀)
    elseif s.t₁ < Inf
        cdf₁,cdf₀     = cdf.([pdist],[τ₁,τ₀])
        ders.ds_dϕ[:] = (ds_dϕ(pdist, τ₁, cdf₁, int_coefs)-ds_dϕ(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.ds_dc    = (ds_dc(pdist, τ₁, cdf₁, int_coefs)-ds_dc(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.ds²_dϕ²  = (ds²_dϕ²(pdist, τ₁, cdf₁, int_coefs)-ds²_dϕ²(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.ds²_dcdϕ = (ds²_dcdϕ(pdist, τ₁, cdf₁, int_coefs)-ds²_dcdϕ(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.ds²_dc²  = (ds²_dc²(pdist, τ₁, cdf₁, int_coefs)-ds²_dc²(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
    else
        cdf₀          = cdf(pdist,τ₀)
        ders.ds_dϕ[:] = -ds_dϕ(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.ds_dc    = -ds_dc(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.ds²_dϕ²  = -ds²_dϕ²(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.ds²_dcdϕ = -ds²_dcdϕ(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.ds²_dc²  = -ds²_dc²(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
    end
end

ds_dϕ(pdist::Distribution, t, x, int_coefs) = clenshaw_asin.(x,int_coefs.ds_dϕ)
ds_dc(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)
ds²_dϕ²(pdist::Distribution, t, x, int_coefs) = clenshaw_asin.(x,int_coefs.ds²_dϕ²)
ds²_dcdϕ(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)*ds_dϕ(pdist, t)
ds²_dc²(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)*ds_dc(pdist, t)

function IntCoefs(pdist::Distribution, degreetype = Val{50}())
    IntCoefs([clenshaw_coefs(pdist, t -> ds_dϕ(pdist,t)[i], degreetype) for i in eachindex(pdist.params)],
    [clenshaw_coefs(pdist, t -> ds²_dϕ²(pdist,t)[i,j], degreetype) for i in eachindex(pdist.params), j in eachindex(pdist.params)])
end

###############################################################
######################GAMMA####################################
###############################################################

type PGamma <: Distribution
    params::Array{Float64}
end

for op in [:pdf, :cdf, :quantile]
    @eval ($op)(pdist::PGamma, t) = ($op)(Gamma(exp(pdist.params[1]),1.), t)
end


ds_dϕ(pdist::PGamma, t) = [log(t)-polygamma(0,pdist.params[1])]
ds_dc(pdist::PGamma, t) = pdist.params[1]-2+t
ds²_dϕ²(pdist::PGamma, t) = reshape([-polygamma(1,pdist.params[1])],1,1)
ds²_dcdϕ(pdist::PGamma, t) = [-1]
ds²_dc²(pdist::PGamma, t) = - t

ds²_dϕ²(pdist::PGamma, t, x, int_coefs) = x*ds²_dϕ²(pdist::PGamma, t)
