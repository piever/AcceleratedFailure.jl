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
        ders.ds_dϕ[:]     = ds_dϕ(pdist, τ₀)
        ders.ds_dc        = ds_dc(pdist, τ₀)
        ders.d²s_dϕ²[:,:] = d²s_dϕ²(pdist, τ₀)
        ders.d²s_dcdϕ[:]  = d²s_dcdϕ(pdist, τ₀)
        ders.d²s_dc²      = d²s_dc²(pdist, τ₀)
    elseif s.t₁ < Inf
        cdf₁,cdf₀          = cdf.([pdist],[τ₁,τ₀])
        ders.ds_dϕ[:]      = (ds_dϕ(pdist, τ₁, cdf₁, int_coefs)-ds_dϕ(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.ds_dc         = (ds_dc(pdist, τ₁, cdf₁, int_coefs)-ds_dc(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)
        ders.d²s_dϕ²[:,:]  = (d²s_dϕ²(pdist, τ₁, cdf₁, int_coefs)-d²s_dϕ²(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)+
                             (ds_dϕds_dϕᵗ(pdist, τ₁, cdf₁, int_coefs)-ds_dϕds_dϕᵗ(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)-
                             ders.ds_dϕ*ders.ds_dϕ'
        ders.d²s_dcdϕ[:]   = (d²s_dcdϕ(pdist, τ₁, cdf₁, int_coefs)-d²s_dcdϕ(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀)-
                              ders.ds_dϕ*ders.ds_dc
        ders.d²s_dc²       = (d²s_dc²(pdist, τ₁, cdf₁, int_coefs)-d²s_dc²(pdist, τ₀, cdf₀, int_coefs))/(cdf₁-cdf₀) -
                              ders.ds_dc*ders.ds_dc
    else
        cdf₀               = cdf(pdist,τ₀)
        ders.ds_dϕ[:]      = -ds_dϕ(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.ds_dc         = -ds_dc(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)
        ders.d²s_dϕ²[:,:]  = -d²s_dϕ²(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)+
                             (-ds_dϕds_dϕᵗ(pdist, τ₀, cdf₀, int_coefs))/(1-cdf₀)-
                             ders.ds_dϕ*ders.ds_dϕ'
        ders.d²s_dcdϕ[:]   = -d²s_dcdϕ(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)-
                              ders.ds_dϕ*ders.ds_dc
        ders.d²s_dc²       = -d²s_dc²(pdist, τ₀, cdf₀, int_coefs)/(1-cdf₀)-
                              ders.ds_dc*ders.ds_dc
    end
end

ds_dϕ(pdist::Distribution, t, x, int_coefs) = clenshaw_asin.(x,int_coefs.ds_dϕ)
ds_dc(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)
d²s_dϕ²(pdist::Distribution, t, x, int_coefs) = clenshaw_asin.(x,int_coefs.d²s_dϕ²)
ds_dϕds_dϕᵗ(pdist::Distribution, t, x, int_coefs) = clenshaw_asin.(x,int_coefs.ds_dϕds_dϕᵗ)
d²s_dcdϕ(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)*ds_dϕ(pdist, t)
d²s_dc²(pdist::Distribution, t, x, int_coefs) = -t*pdf(pdist, t)*ds_dc(pdist, t)


###############################################################
######################GAMMA####################################
###############################################################

type PGamma <: Distribution
    params::Array{Float64}
end

for op in [:pdf, :cdf, :quantile]
    @eval ($op)(pdist::PGamma, t) = ($op)(Gamma(exp(pdist.params[1]),1.), t)
end


ds_dϕ(pdist::PGamma, t)    =  exp(pdist.params[1])*[log(t)-polygamma(0,exp(pdist.params[1]))]
ds_dc(pdist::PGamma, t)    = -exp(pdist.params[1])+t
d²s_dϕ²(pdist::PGamma, t)  = reshape(ds_dϕ(pdist, t)-exp(2*pdist.params[1])*[polygamma(1,exp(pdist.params[1]))],1,1)
d²s_dcdϕ(pdist::PGamma, t) = [-exp(pdist.params[1])]
d²s_dc²(pdist::PGamma, t)  = - t
