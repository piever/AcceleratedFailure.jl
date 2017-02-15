function aft_l(S, X, Xβ, ϕβ, pdist, M, N)
    update_Xβ_dist!(X, Xβ, ϕβ, pdist,  M, N)
    ders = Derivatives(M)
    # compute likelihood
    nll = 0.
    for (i,s) in enumerate(S)
        nll -= compute_loglik!(ders, s, pdist, Xβ[i])
    end
    return nll
end

function aft_h!(grad, hes, S, X, Xβ, ϕβ, pdist, M, N, degreetype)
    update_Xβ_dist!(X, Xβ, ϕβ, pdist, M, N)
    ders = Derivatives(M)
    # Compute auxiliary vector for derivatives.
    int_coefs = IntCoefs(pdist, degreetype)

    # Compute likelihoog, grad and hes

    nll = 0.
    grad[:] = 0.
    hes[:,:] = 0.

    for (i,s) in enumerate(S)
        x = @view X[i,:]
        nll -= compute_loglik!(ders, s, pdist, Xβ[i])
        compute_ders!(ders, s, pdist, Xβ[i], int_coefs)
        subtract_ders!(grad,hes, ders, x, M , N)
    end
    return nll
end

aft{T<:Real}(S::AbstractVector{Event{T}},X::AbstractArray, pdist, degreetype; kwargs...) =
aft(EventWindow.(S),X::AbstractArray, pdist, degreetype; kwargs...)

function aft(S::AbstractVector,X::AbstractArray, pdist, degreetype; kwargs...)
    M = length(pdist.params)
    N = size(X,2)
    Xβ = zeros(size(X,1))
    f1 = (ϕβ) -> aft_l(S, X, Xβ, ϕβ, pdist, M, N)
    h1! = (ϕβ,grad,hes) -> aft_h!(grad, hes, S, X, Xβ, ϕβ, pdist, M, N, degreetype)
    return newton_raphson(f1,h1!, initialize_aft(S, X, pdist, N); kwargs...)
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
    pdist.params[:] = ϕβ[1:length(pdist.params)]
    coefmat = CoefTable(hcat([ϕβ, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], colnames, 4)
    AftModel("Accelerated Failure Time, dist = $pdist;\n", formula, coefmat, M, -neg_ll, -grad, hes, pdist)
end
