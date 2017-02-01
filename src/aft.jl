function aft_l(S, X, ϕβ, dist_type, M, N)
    ϕ = ϕβ[1:M]
    β = ϕβ[(M+1):N]
    Xβ = X*β
    loglik = compute_loglik.(S,[dist_type(ϕ)], Xβ)
    return -sum(loglik)
end

function aft_h!(grad, hes, S, X, ϕβ, dist_type ,M, N, accum_big_ders, degreetype)

    ϕ = ϕβ[1:M]
    β = ϕβ[(M+1):N]
    Xβ = X*β
    zero!(accum_big_ders)
    ders = Derivatives(M)

    # Compute auxiliary vector for derivatives.
    int_coefs = IntCoefs(pdist, degreetype)

    for (i,s) in enumerate(S)
        x = @view X[i,:]
        compute_ders!(ders, s, dist_type(ϕ), Xβ[i], int_coefs)
        add_ders!(accum_big_ders, ders, x, M , N)
    end
    grad[:] = - accum_big_ders.gradlog
    hes[:,:] = - accum_big_ders.heslog
    return - accum_big_ders.val
end


function aft(S::AbstractVector,X::AbstractArray, dist_type, degreetype; kwargs...)
    M = length(pdist.params)
    N = size(X,2)
    accum_big_ders = SmoothLog(0., zeros(M+N), zeros(M+N,M+N))
    #f1 = (β) -> cox_f(S, fs, ls, ξ , X, β, l2_cost, alive, afterΘ) Pensaci
    h1! = (ϕβ,grad,hes) -> aft_h!(grad, hes, S, X, ϕβ, pdist, M, N, accum_big_ders, degreetype)
    return newton_raphson(f1,h1!, vcat(pdist.start,zeros(size(X,2))); kwargs...)
end

function aft(formula::Formula, data::DataFrame, dist_type, degreetype = Val{50}(); kwargs...)
    sorted_data = sort(data, cols = formula.lhs)
    M = DataFrames.ModelFrame(formula,sorted_data)
    S = convert(Array, M.df[:,1])
    model_matrix = DataFrames.ModelMatrix(M)
    X = model_matrix.m
    β, neg_ll,grad, hes =  aft(S, X, dist_type, degreetype; kwargs...)
    colnames = coefnames(M)
    se = sqrt.(diag(pinv(hes)))
    z_score = β./se
    pvalues = 2*cdf(Normal(),-abs.(z_score))
    coefmat = CoefTable(hcat([β, se, z_score, pvalues]...),
    ["Estimate", "Std.Error", "z value", "Pr(>|z|)"], colnames, 4)
    EventHistoryModel("Accelerated Failure Time, dist = $pdist;\n", formula, coefmat, M, -neg_ll, -grad, hes)
end
