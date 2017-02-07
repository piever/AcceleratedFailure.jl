function aft_l(S, X, ϕβ, pdist, M, N)
    ϕ = ϕβ[1:M]
    β = ϕβ[M+(1:N)]
    Xβ = X*β
    # Update pdist
    pdist.params[:] = ϕ
    # compute likelihood
    loglik = compute_loglik.(S,[pdist], Xβ)
    return -sum(loglik)
end

function aft_h!(grad, hes, S, X, ϕβ, pdist ,M, N, accum_big_ders, degreetype)

    ϕ = ϕβ[1:M]
    β = ϕβ[M+(1:N)]
    Xβ = X*β
    zero!(accum_big_ders)
    ders = Derivatives(M)

    # Update pdist
    pdist.params[:] = ϕ
    # Compute auxiliary vector for derivatives.
    int_coefs = IntCoefs(pdist, degreetype)

    for (i,s) in enumerate(S)
        x = @view X[i,:]
        compute_ders!(ders, s, pdist, Xβ[i], int_coefs)
        add_ders!(accum_big_ders, ders, x, M , N)
    end
    grad[:] = - accum_big_ders.gradlog
    hes[:,:] = - accum_big_ders.heslog
    return - accum_big_ders.val
end

aft{T<:Real}(S::AbstractVector{Event{T}},X::AbstractArray, pdist, degreetype; kwargs...) =
aft(EventWindow.(S),X::AbstractArray, pdist, degreetype; kwargs...)

function aft(S::AbstractVector{EventWindow},X::AbstractArray, pdist, degreetype; kwargs...)
    St = [(s.t₀+s.t₁)/2 for s in S]
    valid = isfinite(St)
    res = glm(X[valid,:], St[valid], Gamma(), LogLink())
    T = 1/GLM.dispersion(res, true)
    M = length(pdist.params)
    N = size(X,2)
    accum_big_ders = SmoothLog(0., zeros(M+N), zeros(M+N,M+N))
    f1 = (ϕβ) -> aft_l(S, X, ϕβ, pdist, M, N)
    h1! = (ϕβ,grad,hes) -> aft_h!(grad, hes, S, X, ϕβ, pdist, M, N, accum_big_ders, degreetype)
    return newton_raphson(f1,h1!, vcat(log(T),coef(res)); kwargs...)
end

function aft(formula::Formula, data::DataFrame, pdist, degreetype = Val{50}(); kwargs...)
    M = DataFrames.ModelFrame(formula,data)
    S = convert(Array, M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m
    ϕβ, neg_ll,grad, hes =  aft(S, X, pdist, degreetype; kwargs...)
    colnames = vcat(string.(["params"],1:length(pdist.params)),coefnames(M))
    se = sqrt.(diag(pinv(hes)))
    z_score = ϕβ./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([ϕβ, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], colnames, 4)
    SurvivalModel("Accelerated Failure Time, dist = $pdist;\n", formula, coefmat, M, -neg_ll, -grad, hes)
end
