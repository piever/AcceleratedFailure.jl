function compute_loglik(s, dist, coef)
    if s.t₁ == s.t₀
        return log(pdf(dist,s.t₀*exp(-coef)))-coef
    elseif s.t₁ < Inf
        return log(cdf(dist,s.t₁*exp(-coef))-cdf(dist,s.t₀*exp(-coef)))
    else
        return log(1-cdf(dist,s.t₀*exp(-coef)))
    end
end

function compute_ders!(ders, s, pdist::Distribution, c, int_coefs)
    M = length(pdist.params)
    if s.t₁ == s.t₀
        τ₀ = s.t₀*exp(-c)
        ders.loglik       = log(pdf(pdist,τ₀))-c
        for i in 1:M
            ders.ds_dϕ[i] = ds_dϕ(pdist, τ₀, int_coefs,i)
        end
        ders.ds_dc        = ds_dc(pdist, τ₀, int_coefs)
        for j in 1:M, i in 1:M
            ders.d²s_dϕ²[i,j] = d²s_dϕ²(pdist, τ₀,int_coefs, i, j)
        end
        for i in 1:M
            ders.d²s_dcdϕ[i]  = d²s_dcdϕ(pdist, τ₀, int_coefs, i)
        end
        ders.d²s_dc²      = d²s_dc²(pdist, τ₀, int_coefs)
    elseif s.t₁ < Inf
        τ₁,τ₀              = [s.t₁,s.t₀].*exp(-c)
        cdf₁,cdf₀          = cdf.([pdist],[τ₁,τ₀])
        p₀                 = τ₀*pdf(pdist, τ₀)
        p₁                 = τ₁*pdf(pdist, τ₁)
        ders.loglik        = log(cdf₁-cdf₀)
        for i in 1:M
            ders.ds_dϕ[i]      = (ds_dϕ_int(pdist, τ₁, cdf₁, int_coefs,i)-ds_dϕ_int(pdist, τ₀, cdf₀, int_coefs,i))/(cdf₁-cdf₀)
        end
        ders.ds_dc         = -(p₁-p₀)/(cdf₁-cdf₀)
        for j in 1:M, i in 1:M
            ders.d²s_dϕ²[i,j]  = (d²s_dϕ²_int(pdist, τ₁, cdf₁, int_coefs,i,j)-d²s_dϕ²_int(pdist, τ₀, cdf₀, int_coefs,i,j))/(cdf₁-cdf₀)-
                                 ders.ds_dϕ[i]*ders.ds_dϕ[j]
        end
        for i in 1:M
            ders.d²s_dcdϕ[i]   = -(p₁*ds_dϕ(pdist, τ₁, int_coefs, i)-p₀*ds_dϕ(pdist, τ₀, int_coefs, i))/(cdf₁-cdf₀)-
                                  ders.ds_dϕ[i]*ders.ds_dc
        end
        ders.d²s_dc²       = -(p₁*ds_dc(pdist, τ₁, int_coefs)-p₀*ds_dc(pdist, τ₀, int_coefs))/(cdf₁-cdf₀) -
                              ders.ds_dc*ders.ds_dc
    else
        τ₀ = s.t₀*exp(-c)
        cdf₀               = cdf(pdist,τ₀)
        p₀                 = τ₀*pdf(pdist, τ₀)
        ders.loglik        = log(1-cdf₀)
        for i in 1:M
            ders.ds_dϕ[i]      = -ds_dϕ_int(pdist, τ₀, cdf₀, int_coefs,i)/(1-cdf₀)
        end
        ders.ds_dc         = p₀/(1-cdf₀)
        for j in 1:M, i in 1:M
            ders.d²s_dϕ²[i,j]  = -d²s_dϕ²_int(pdist, τ₀, cdf₀, int_coefs,i,j)/(1-cdf₀)-
                                ders.ds_dϕ[i]*ders.ds_dϕ[j]
        end
        for i in 1:M
            ders.d²s_dcdϕ[i]   = p₀*ds_dϕ(pdist, τ₀, int_coefs, i)/(1-cdf₀)-
                                  ders.ds_dϕ[i]*ders.ds_dc
        end
        ders.d²s_dc²       = p₀*ds_dc(pdist, τ₀, int_coefs)/(1-cdf₀)-
                              ders.ds_dc*ders.ds_dc
    end
end

ds_dϕ_int(pdist::Distribution, t, x, int_coefs, i) = clenshaw_asin(x,int_coefs.g[i])
d²s_dϕ²_int(pdist::Distribution, t, x, int_coefs, i, j) = clenshaw_asin(x,int_coefs.hggt[i,j])

auxvec(pdist::Distribution) = Array{Float64,1}(0)

initialize_aft(S::AbstractVector,X::AbstractArray, pdist::Distribution, N) = vcat(pdist.params,zeros(N))
###############################################################
######################GAMMA####################################
###############################################################

type PGamma <: Distribution
    params::Array{Float64}
end

PGamma() = PGamma([1.])

for op in [:pdf, :cdf, :quantile]
    @eval ($op)(pdist::PGamma, t) = ($op)(Gamma(exp(pdist.params[1]),exp(-pdist.params[1])), t)
end

rand(pdist::PGamma) = rand(Gamma(exp(pdist.params[1]),exp(-pdist.params[1])))

ds_dϕ(pdist::PGamma, t, int_coefs, i)    =  exp(pdist.params[1])*(log(t)-t) + int_coefs.aux[1]
ds_dc(pdist::PGamma, t, int_coefs)    = -exp(pdist.params[1])*(1-t)
d²s_dϕ²(pdist::PGamma, t, int_coefs,i, j)  = ds_dϕ(pdist, t, int_coefs,i)+ int_coefs.aux[2]
d²s_dcdϕ(pdist::PGamma, t, int_coefs, i) = -exp(pdist.params[1])*(1-t)
d²s_dc²(pdist::PGamma, t, int_coefs)  = - exp(pdist.params[1])*t

auxvec(pdist::PGamma) = [exp(pdist.params[1])*(-polygamma(0,exp(pdist.params[1]))+1+pdist.params[1]),
-exp(2*pdist.params[1])*polygamma(1,exp(pdist.params[1]))+exp(pdist.params[1])]

function initialize_aft(S::AbstractVector,X::AbstractArray, pdist::PGamma, N)
    St = [(s.t₀+s.t₁)/2 for s in S]
    valid = isfinite(St)
    res = glm(X[valid,:], St[valid], Gamma(), LogLink())
    T = 1/GLM.dispersion(res, true)
    return vcat(log(T),coef(res))
end
