surv(f::Function, pts::LinSpace) = surv(f.(pts),pts)
function surv(v,pts::LinSpace)
    binsize = step(pts)
    w = cumsum(v)*binsize
    w[:] = w[end]-w+binsize*v/2
    return w
end

interpolate(x, f::Function, pts::LinSpace) = interpolate(x, f.(pts),pts)
function interpolate(x,v,pts::LinSpace)
    if x <= pts.start
        return v[1]
    elseif x >= pts.stop
        return v[end]
    end
    ratio = (x-pts.start)/step(pts)
    rem = ratio-floor(ratio)
    I1 = floor(Int,ratio)+1
    value = v[I1]*(1.-rem)+v[I1+1]*rem
    return value
end


# τ = t*exp(-dot(β,x)) where t is measured time and τ is baseline time
#inputs =  grad and hessian of F(ϕ,t), length of ϕ, length of β,β,τ,x
# output = grad and hessian of F(ϕ,t*exp(-dot(β,x))) wrt ϕ and β
function aft_gradhes(dFdϕτ, d2Fdϕτ2,τ,x,M,N,β)
    dτdβ = -τ*x
    d2τdβ2 = τ*x*x'
    grad = zeros(M+N)
    grad[1:M] = dFdϕτ[1:M]
    grad[M+(1:N)] = dFdϕτ[M+1]*dτdβ
    hes = zeros(M+N,M+N)
    hes[1:M,1:M] = d2Fdϕτ2[1:M,1:M]
    hes[1:M, M+(1:N)] = d2Fdϕτ2[1:M,M+1]*dτdβ'
    hes[M+(1:N),1:M] = hes[1:M, M+(1:N)]'
    hes[M+(1:N),M+(1:N)] = dFdϕτ[M+1]*d2τdβ2+d2Fdϕτ2[M+1,M+1]*dτdβ*dτdβ'
    return grad,hes
end




# careful! eff_t0 could go to 0!!!
function aft_h!(grad, hes, P, P1, P2, ϕ, ts, X,M,N, β)
    #eff_t0 = max.(times.*exp(-X*β),1e-3)
    Xβ = X*β
    τs = ts.*exp(-Xβ)
    ll = sum(log.(pdf(Gamma(T,1.),eff_t0[cs]))-Xβ)
    grad[:] = 0.
    hess[:,:] = 0.
    for i = 1:length(ts)
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




# # careful! eff_t0 could go to 0!!!
# function get_likelihood(t0, cs,t1, X, T, β)
#     #eff_t0 = max.(times.*exp(-X*β),1e-3)
#     Xβ = X*β
#     eff_t0 = t0.*exp(-Xβ)
#     eff_t1 = t1.*exp(-Xβ)
#     ll = sum(log.(pdf(Gamma(T,1.),eff_t0[cs]))-Xβ)+
#         sum(log.(cdf(Gamma(T,1.),eff_t1[cs])-cdf(Gamma(T,1.),eff_t1[cs])))
#     # get gradient and hessian
#     # bin x axis (to compute gradient and hessian for cdf)
#     ϵ = 1e-4
#     nbins = 500
#     disc_axis = linspace(0.,quantile(Gamma(T,1.),1.-ϵ),nbins)
#     gammas = pdf(Gamma(T,1), disc_axis)
#     # compute first two derivatives of 1-cdf as fct of threshold
#     dsurvdT = surv(dlgammadT.(T,disc_axis).*gammas,disc_axis)
#     dsurvdT2 = surv(dlgammadT2.(T,disc_axis).*gammas,disc_axis)
#     # initialize gradient and hessian
#     grad = zeros(1,1+size(X,2))
#     hes = zeros(1+size(X,2),1+size(X,2))
#     for i in 1:length(eff_t0)
#         addgrad!(grad, T, eff_t0[i],cs[i], X[i:i,:],dsurvdT,disc_axis)
#         addhes!(hes, T, eff_t0[i], cs[i], X[i:i,:],dsurvdT,dsurvdT2,disc_axis)
#     end
#     return -ll, -grad, -hes
# end
