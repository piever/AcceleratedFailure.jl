function cox_f(S, fs, ls, ξ , X, β, λ, alive, afterΘ) # preprocessed already
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!
    Xβ = X*β
    Θ = exp.(Xβ)
    after!(afterΘ, Θ)

    y = 0.
    #compute loglikelihood, score, fischer_info
    for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterΘ[fs[i]]-ρ*(afterΘ[fs[i]]-afterΘ[ls[i]+1])
            y -= Xβ[j] -log(ϕ)
        end
    end
    y += λ*dot(β',β)
    return y
end

# preprocessed already:
# fs = index first deaths, ls = index last deaths,
# X is covariates, ξ is covariate covariate transpose
function cox_h!(grad,hes, S, fs, ls, ξ , X, β, λ, alive, afterΘ, afterXΘ, afterξΘ)
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!

    Xβ = X*β
    Θ = exp.(Xβ)
    after!(afterΘ,Θ)
    after!(afterXΘ, X.*Θ)
    after!(afterξΘ,ξ.*Θ)

    y = 0.
    grad[:] = 0.
    hes[:] = 0.

    # preallocate
    Z = zeros(size(X,2))
    Ξ = zeros(size(X,2),size(X,2))

    #compute loglikelihood, score, fischer_info
    for i in 1:length(fs)
        for j in (fs[i]):(ls[i])
            ρ = (alive[j]-alive[fs[i]])/(alive[fs[i]]-alive[ls[i]+1])
            ϕ = afterΘ[fs[i]]-ρ*(afterΘ[fs[i]]-afterΘ[ls[i]+1])
            Z[:] = afterXΘ[fs[i],:]-ρ*(afterXΘ[fs[i],:]-afterXΘ[ls[i]+1,:])
            Ξ[:,:] = afterξΘ[fs[i],:,:]-ρ*(afterξΘ[fs[i],:,:]-afterξΘ[ls[i]+1,:,:])
            y -= Xβ[j] -log(ϕ)
            grad[:] += - X[j,:]+Z/ϕ
            hes[:, :] += Ξ/ϕ - Z*Z'/ϕ^2
        end
    end
    y += λ*dot(β',β)
    grad[:] +=  2*λ*β
    hes[:,:] +=  2*λ*eye(size(X,2))
    return y
end

function coxph(S::AbstractVector,X::AbstractArray; l2_cost = 0., kwargs...)
    ξ = zeros(size(X,1),size(X,2),size(X,2))
    for i in 1:size(X,1)
        ξ[i,:,:] = X[i,:]*X[i,:]'
    end

    # compute first and last!

    fs = find(firsts(S))
    ls = find(lasts(S))
    # do optimization

    alive = after(ones(Int64, length(S)))
    afterΘ = init_after(zeros(size(X,1)))
    afterXΘ = init_after(X.*zeros(size(X,1)))
    afterξΘ = init_after(ξ.*zeros(size(X,1)))

    f1 = (β) -> cox_f(S, fs, ls, ξ , X, β, l2_cost, alive, afterΘ)
    h1! = (β,grad,hes) -> cox_h!(grad,hes, S, fs, ls, ξ , X, β, l2_cost,alive, afterΘ, afterXΘ, afterξΘ)
    return newton_raphson(f1,h1!, zeros(size(X,2)); kwargs...)
end

function coxph(formula::Formula, data::DataFrame; l2_cost = 0., kwargs...)
    sorted_data = sort(data, cols = formula.lhs)
    #sort!(sorted_data, cols = formula.lhs)
    M = DataFrames.ModelFrame(formula,sorted_data)
    S = convert(Array, M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m[:,2:size(model_matrix.m,2)]
    β, neg_ll,grad, hes =  coxph(S, X; l2_cost = l2_cost, kwargs...)
    colnames = coefnames(M)[2:end]
    se = sqrt.(diag(pinv(hes)))
    z_score = β./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], colnames, 4)
    EventHistoryModel("Cox", formula, coefmat, M, -neg_ll, hes)
end
