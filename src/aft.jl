function get_likelihood(S, X, ϕ, β, pdist)
    Xβ = X*β
    coef = exp(-Xβ)
    loglik = sum(ll.(S,[pdist],[ϕ], coef))
    return loglik
end
