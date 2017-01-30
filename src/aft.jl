function aft_l(S, X, ϕ, β, pdist)
    Xβ = X*β
    coef = exp(-Xβ)
    loglik = sum(ll.(S,[pdist],[ϕ], coef))
    return loglik
end


# careful! eff_t0 could go to 0!!!
function aft_h!(grad, hes, S, ϕ, β, pdist, X,N)
    #eff_t0 = max.(times.*exp(-X*β),1e-3)
    Xβ = X*β
    coef = exp(-Xβ)
    τs = ts.*coef
    ll = sum(ll.(S,[pdist],[ϕ], coef))








    grad[:] = 0.
    hess[:,:] = 0.
    for i = 1:length(S)

        dFdϕτ = P1(ϕ,τs[i])
        dF2dϕτ2 = P2(ϕ,τs[i])
        gradt,hest = aft_gradhes(dFdϕτ, d2Fdϕτ2,τs[i],X[i,:],M,N,β)
        gradt -= X[i,:]
        grad += gradt
        hes += hest
    end
    grad[:] = -grad
    hes[:,:] = -hes
    return -ll
end
