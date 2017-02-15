function update_Xβ_dist!(X, Xβ, ϕβ, pdist, M, N)
    ϕ = ϕβ[1:M]
    β = ϕβ[M+1:end]
    A_mul_B!(Xβ, X, β)
    pdist.params[:] = ϕ
    return
end

function compute_loglik!(ders, s, pdist::Distribution, c)
    if s.t₁ == s.t₀
        ders.τs[1] = s.t₀*exp(-c)
        return log(pdf(pdist, ders.τs[1]))-c
    else
        fin = s.t₁ < Inf
        for k in 1:(fin ? 2 : 1)
            ders.τs[k]                 = ((k ==1) ? s.t₀ : s.t₁)*exp(-c)
            ders.cdfs[k]               = cdf(pdist, ders.τs[k])
        end
        ders.Δcdf[1]               = (fin ? ders.cdfs[2] : 1.) - ders.cdfs[1]
        return log(ders.Δcdf[1])
    end
end

# Add other fields to ders object for in-place operations
# Clean this bit
function compute_ders!(ders, s, pdist::Distribution, c, int_coefs)
    M = length(pdist.params)
    if s.t₁ == s.t₀
        d²l!(ders.heslog, ders.gradlog, pdist, ders.τs[1], int_coefs)
    else
        fin = s.t₁ < Inf
        for k in 1:(fin ? 2 : 1)
            ders.ps[k]                 = ders.τs[k]*pdf(pdist, ders.τs[k])
            d²l_int!(ders.gradlogint[k], ders.heslogint[k], pdist, ders.τs[k], ders.cdfs[k], ders.ps[k], int_coefs)
        end
        for i in eachindex(ders.gradlog)
            ders.gradlog[i] = ((fin ?  ders.gradlogint[2][i] : 0) - ders.gradlogint[1][i])/ders.Δcdf[1]
        end
        for j in 1:M+1, i in 1:M+1
            ders.heslog[i,j] = ((fin ? ders.heslogint[2][i,j] : 0) - ders.heslogint[1][i,j])/ders.Δcdf[1] -
                                ders.gradlog[i]*ders.gradlog[j]
        end
    end
end

function subtract_ders!(grad, hes, ders::Derivatives, x, M , N)
    for i = 1:(M+N)
        grad[i] -= ders.gradlog[min(i,M+1)]*x[max(i-M,1)]
    end
    for j = 1:(M+N), i = 1:(M+N)
        hes[i,j] -= ders.heslog[min(i,M+1),min(j,M+1)]*x[max(i-M,1)]*x[max(j-M,1)]
    end
end

function d²l_int!(gl, hl, pdist::Distribution, t, x, p, int_coefs)
    for j1 = 1:length(pdist.params)
        gl[j1] = clenshaw_asin(x,int_coefs.g[j1])
    end
    gl[end] = -p
    last_line = @view hl[end, :]
    dl!(last_line, pdist, t, int_coefs)
    scale!(last_line, -p)
    for j2 in 1:length(pdist.params), j1 in j2:length(pdist.params)
        hl[j1,j2] = clenshaw_asin(x,int_coefs.hggt[j1,j2])
    end
    for j2 in 1:length(pdist.params)+1, j1 in 1:j2-1
        hl[j1,j2] = hl[j2,j1]
    end
end


auxvec(pdist::Distribution) = Array{Float64,1}(0)

initialize_aft(S::AbstractVector,X::AbstractArray, pdist::Distribution, N) = vcat(pdist.params,zeros(N))
