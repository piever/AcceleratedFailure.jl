function update_Xβ_dist!(rr, ϕβ)
    ϕ = ϕβ[1:rr.M]
    β = ϕβ[rr.M+1:end]
    A_mul_B!(rr.Xβ, rr.X, β)
    rr.Θ .= exp.(-rr.Xβ)
    rr.dist.params[:] = ϕ
    return
end

function compute_loglik!(rr)
    for i in rr.I1
        rr.τs[i,1] = rr.W[i,1]*rr.Θ[i]
        rr.loglik[i] = logpdf(rr.dist, rr.τs[i,1])-rr.Xβ[i]
    end
    for i in rr.I2
        for k in 1:(rr.W[i,2] < Inf ? 2 : 1)
            rr.τs[i, k]    = rr.W[i,k]*rr.Θ[i]
            rr.cdfs[i, k]  = cdf(rr.dist, rr.τs[i,k])
        end
        rr.Δcdf[i] = rr.cdfs[i,2] - rr.cdfs[i,1]
        rr.loglik[i] = log(rr.Δcdf[i])
    end
    return sum(rr.loglik)
end

function diff_ders!(gl, hl, gli, hli, fin, Δcdf)
    @inbounds for j1 in eachindex(gl)
        gl[j1] = (fin*gli[2][j1]  - gli[1][j1])/Δcdf
    end
    @inbounds for j2 in eachindex(gl), j1 in eachindex(gl)
        hl[j1,j2] = (fin*hli[2][j1,j2]- hli[1][j1,j2])/Δcdf - gl[j1]*gl[j2]
    end
end

function subtract_ders!(grad, hes, rr, int_coefs)
    grad[:] = 0.
    hes[:,:] = 0.
    gl, hl = zeros(rr.M+1), zeros(rr.M+1,rr.M+1)
    gli, hli = [zeros(rr.M+1) for i in 1:2], [zeros(rr.M+1,rr.M+1) for i in 1:2]
    for i in rr.I1
        d²l!(gl, hl, rr.dist, rr.τs[i,1], int_coefs)
        subtract_ders!(grad, hes, gl, hl, @view(rr.X[i,:]), rr.M, rr.N)
    end
    for i in rr.I2
        fin = rr.W[i,2] < Inf
        @inbounds for k in 1:(fin ? 2 : 1)
            d²l_int!(gli[k], hli[k], rr.dist, rr.τs[i,k], rr.cdfs[i,k], int_coefs)
        end
        diff_ders!(gl, hl, gli, hli, fin, rr.Δcdf[i])
        subtract_ders!(grad, hes, gl, hl, @view(rr.X[i,:]), rr.M, rr.N)
    end
end

function subtract_ders!(grad, hes, gl, hl, x, M , N)
    @inbounds for j1 = 1:(M+N)
        grad[j1] -= gl[min(j1,M+1)]*x[max(j1-M,1)]
    end
    @inbounds for j2 = 1:(M+N), j1 = 1:(M+N)
        hes[j1,j2] -= hl[min(j1,M+1),min(j2,M+1)]*x[max(j1-M,1)]*x[max(j2-M,1)]
    end
end

function d²l_int!(gl, hl, pdist::Distribution, t, x, int_coefs)
    @inbounds for j1 = 1:length(pdist.params)
        gl[j1] = clenshaw_asin(x,int_coefs.g[j1])
    end
    gl[end] = -t*pdf(pdist, t)
    last_line = @view hl[end, :]
    dl!(last_line, pdist, t, int_coefs)
    scale!(last_line, gl[end])
    @inbounds for j2 in 1:length(pdist.params), j1 in j2:length(pdist.params)
        hl[j1,j2] = clenshaw_asin(x,int_coefs.hggt[j1,j2])
    end
    @inbounds for j2 in 1:length(pdist.params)+1, j1 in 1:j2-1
        hl[j1,j2] = hl[j2,j1]
    end
end


auxvec(pdist::Distribution) = Array{Float64,1}(0)

initialize_aft(S::AbstractVector,X::AbstractArray, pdist::Distribution, N) =
vcat(pdist.params,zeros(N))
