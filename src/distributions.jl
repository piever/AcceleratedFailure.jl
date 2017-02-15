###############################################################
######################GAMMA####################################
###############################################################

immutable PGamma <: Distribution
    params::Array{Float64}
end

PGamma() = PGamma([1.])

for op in [:pdf, :cdf, :quantile]
    @eval ($op)(pdist::PGamma, t) = ($op)(Gamma(exp(pdist.params[1]),exp(-pdist.params[1])), t)
end

rand(pdist::PGamma) = rand(Gamma(exp(pdist.params[1]),exp(-pdist.params[1])))

ds(pdist::PGamma, t, int_coefs, i)    =  exp(pdist.params[1])*(log(t)-t) + int_coefs.aux[1]
d²s(pdist::PGamma, t, int_coefs,i, j)  = ds(pdist, t, int_coefs,i)+ int_coefs.aux[2]

function dl!(gl, pdist::PGamma, t, int_coefs)
    gl[2] = exp(pdist.params[1])
    gl[1] = muladd(gl[2],(log(t)-t), int_coefs.aux[1])
    gl[2] *= (t-1)
end

function d²l!(hl, gl, pdist::PGamma, t, int_coefs)
    hl[2,2] = gl[2] = exp(pdist.params[1])
    gl[1] = muladd(gl[2],(log(t)-t), int_coefs.aux[1])
    gl[2] *= (t-1)
    hl[1,1] = gl[1] + int_coefs.aux[2]
    hl[1,2] = hl[2,1] = gl[2]
    hl[2,2] *= -t
end


auxvec(pdist::PGamma) = [exp(pdist.params[1])*(-polygamma(0,exp(pdist.params[1]))+1+pdist.params[1]),
-exp(2*pdist.params[1])*polygamma(1,exp(pdist.params[1]))+exp(pdist.params[1])]

function initialize_aft(S::AbstractVector,X::AbstractArray, pdist::PGamma, N)
    St = [(s.t₀+s.t₁)/2 for s in S]
    valid = isfinite(St)
    res = glm(X[valid,:], St[valid], Gamma(), LogLink())
    T = 1/GLM.dispersion(res, true)
    return vcat(log(T),coef(res))
end
