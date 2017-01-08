function after{N,T}(v::AbstractArray{T,N})
    cv = cumsum(v,1)
    newsize = collect(size(cv))
    newsize[1] += 1
    afterv = zeros(eltype(cv),newsize...)
    if N==1
        afterv[1:end-1] = cv[end] - cv + v
    elseif N==2
        afterv[1:end-1,:] = (cv[end:end,:] .- cv) + v
    elseif N==3
        afterv[1:end-1,:,:] = (cv[end:end,:,:] .- cv) + v
    else
        error("$v can be at most 3-dimensional!")
    end
    return afterv
end

function cox_f(times, censored, f, l, ξ , X, β, λ) # preprocessed already
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!
    Xβ = X*β
    Θ = exp.(Xβ)
    afterΘ = after(Θ)
    alive = after(ones(Int64, length(times)))

    y = 0.
    #compute loglikelihood, score, fischer_info
    for i in 1:length(f)
        for j in (f[i]):(l[i])
            ρ = (alive[j]-alive[f[i]])/(alive[f[i]]-alive[l[i]+1])
            ϕ = afterΘ[f[i]]-ρ*(afterΘ[f[i]]-afterΘ[l[i]+1])
            y -= Xβ[j] -log(ϕ)
        end
    end
    return y
end

# preprocessed already:
# f = index first deaths, l = index last deaths,
# X is covariates, ξ is covariate covariate transpose
function cox_h!(grad,hes, times, censored, f, l, ξ , X, β, λ)
    #compute relevant quantities for loglikelihood, score, fischer_info
    ## trick= usa scale invariance e semplifica tutto!!!

    Xβ = X*β
    Θ = exp.(Xβ)
    afterΘ = after(Θ)
    alive = after(ones(Int64, length(times)))
    afterXΘ = after(X.*Θ)
    afterξΘ = after(ξ.*Θ)

    y = 0.
    grad[:] = 0.
    hes[:] = 0.
    #compute loglikelihood, score, fischer_info
    for i in 1:length(f)
        for j in (f[i]):(l[i])
            ρ = (alive[j]-alive[f[i]])/(alive[f[i]]-alive[l[i]+1])
            ϕ = afterΘ[f[i]]-ρ*(afterΘ[f[i]]-afterΘ[l[i]+1])
            Z = afterXΘ[f[i],:]-ρ*(afterXΘ[f[i],:]-afterXΘ[l[i]+1,:])
            Ξ = afterξΘ[f[i],:,:]-ρ*(afterξΘ[f[i],:,:]-afterξΘ[l[i]+1,:,:])
            y -= Xβ[j] -log(ϕ)
            grad[:] = grad - X[j,:]+Z/ϕ
            hes[:] = hes + Ξ/ϕ - Z*Z'/ϕ^2
        end
    end
    grad[:] = grad + λ
    hes[:] = hes + λ*eye(size(X,2))
    return y
end

function coxph(S::AbstractVector,X::AbstractArray; penalized = 0.)
    times = [a.time for a in S]
    censored = [a.censored for a in S]
    ξ = zeros(size(X,1),size(X,2),size(X,2))
    for i in 1:size(X,1)
        ξ[i,:,:] = X[i,:]*X[i,:]'
    end

    # compute first and last!
    firsts = !censored & [t==1 || S[t] > S[t-1] for t = 1:length(S)]
    lasts = !censored & [t==length(S) || S[t+1] > S[t] for t = 1:length(S)]

    f = find(firsts)
    l = find(lasts)
    # do optimization

    f1 = (β) -> cox_f(times, censored, f, l, ξ , X, β, penalized)
    h1! = (β,grad,hes) -> cox_h!(grad,hes,times, censored, f, l, ξ , X, β, penalized)
    return newtonraphson(f1,h1!, zeros(size(X,2)))
end

function coxph(formula::Formula, data::DataFrame)
    sorted_data = deepcopy(data)
    sort!(sorted_data, cols = formula.lhs)
    M = DataFrames.ModelFrame(formula,sorted_data)
    S = M.df[:,1]
    X = DataFrames.ModelMatrix(M)
    X = X.m[:,collect(2:size(X.m,2))]
    return coxph(S, X)
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
