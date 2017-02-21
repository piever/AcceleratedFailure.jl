function aft_l(ϕβ, rr)
    update_Xβ_dist!(rr, ϕβ)
    # compute likelihood
    nll = -compute_loglik!(rr)
    return nll
end

function aft_h!(grad, hes, ϕβ, rr, degreetype)
    update_Xβ_dist!(rr, ϕβ)
    # Compute auxiliary vector for derivatives.
    int_coefs = IntCoefs(rr.dist, degreetype)

    # Compute likelihoog, grad and hes
    nll = -compute_loglik!(rr)
    subtract_ders!(grad,hes, rr, int_coefs)
    return nll
end

aft{T<:Real}(S::AbstractVector{Event{T}},X::AbstractArray, pdist, degreetype; kwargs...) =
aft(EventWindow.(S),X::AbstractArray, pdist, degreetype; kwargs...)

function aft(S::AbstractVector,X::AbstractArray, pdist, degreetype; kwargs...)
    rr = AftResp(S, X, pdist)
    f1 = (ϕβ) -> aft_l(ϕβ, rr)
    h1! = (ϕβ,grad,hes) -> aft_h!(grad, hes, ϕβ, rr, degreetype)
    return newton_raphson(f1,h1!, initialize_aft(S, X, pdist, size(X,2)); kwargs...)
end

function aft(formula::Formula, data::DataFrame, pdist, degreetype = Val{50}(); kwargs...)
    M = DataFrames.ModelFrame(formula,data)
    S = convert(Array, M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m
    ϕβ, neg_ll,grad, hes =  aft(S, X, pdist, degreetype; kwargs...)
    rownms = vcat(string.(["params"],1:length(pdist.params)),coefnames(M))
    se = sqrt.(diag(pinv(hes)))
    z_score = ϕβ./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    pdist.params[:] = ϕβ[1:length(pdist.params)]
    coefmat = CoefTable(hcat([ϕβ, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], rownms, 4)
    AftModel("Accelerated Failure Time, dist = $pdist;\n", formula,
             coefmat, M, -neg_ll, -grad, hes, pdist)
end
