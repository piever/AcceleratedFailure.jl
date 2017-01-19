function cox_f(S, fs, ls, ξ , X, β, λ) # preprocessed already
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!
    Xβ = X*β
    Θ = exp.(Xβ)
    afterΘ = after(Θ)
    alive = after(ones(Int64, length(S)))

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
# f = index first deaths, l = index last deaths,
# X is covariates, ξ is covariate covariate transpose
function cox_h!(grad,hes, S, fs, ls, ξ , X, β, λ)
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!

    Xβ = X*β
    Θ = exp.(Xβ)
    afterΘ = after(Θ)
    alive = after(ones(Int64, length(S)))
    afterXΘ = after(X.*Θ)
    afterξΘ = after(ξ.*Θ)

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

    f1 = (β) -> cox_f(S, fs, ls, ξ , X, β, l2_cost)
    h1! = (β,grad,hes) -> cox_h!(grad,hes, S, fs, ls, ξ , X, β, l2_cost)
    return newtonraphson(f1,h1!, zeros(size(X,2)); kwargs...)
end

function coxph(formula::Formula, data::DataFrame; l2_cost = 0., kwargs...)
    sorted_data = sort(data, cols = formula.lhs)
    #sort!(sorted_data, cols = formula.lhs)
    M = DataFrames.ModelFrame(formula,sorted_data)
    S = collect(M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m[:,2:size(model_matrix.m,2)]
    β, neg_ll,grad, hes =  coxph(S, X; l2_cost = l2_cost, kwargs...)
    colnames = coefnames(M)[2:end]
    se = sqrt.(diag(pinv(hes)))
    z_score = β./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], colnames, 4)
    EventHistoryModel("Cox", formula, coefmat, M)
end

# function phreg(formula::Formula, data::DataFrame; id=[], opt...)
#     ## fm = DataFrames.Formula(formula.args[2],formula.args[3])
#     M = DataFrames.ModelFrame(formula,data)
#     S = M.df[:,1] ## convert(Vector,M.df[:,1])
#     time = EventHistory.Time(S)
#     status = EventHistory.Status(S)
#     entry = try EventHistory.Entry(S) catch [] end
#     X = DataFrames.ModelMatrix(M)
#     X = X.m[:,collect(2:size(X.m,2))]
#     res = EventHistory.phreg(X, time, status, entry; opt...)
#     cnames = setdiff(coefnames(M),["(Intercept)"])
#     res.coefmat.rownms=cnames
#     res.formula = formula
#     res.eventtype = typeof(S[1])
#     res
# end
