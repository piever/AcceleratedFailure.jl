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
# careful! eff_t0 could go to 0!!!
function get_likelihood(t0, cs,t1, X, T, β)
    #eff_t0 = max.(times.*exp(-X*β),1e-3)
    Xβ = X*β
    eff_t0 = t0.*exp(-Xβ)
    eff_t1 = t1.*exp(-Xβ)
    ll = sum(log.(pdf(Gamma(T,1.),eff_t0[cs]))-Xβ)+
        sum(log.(cdf(Gamma(T,1.),eff_t1[cs])-cdf(Gamma(T,1.),eff_t1[cs])))
    # get gradient and hessian
    # bin x axis (to compute gradient and hessian for cdf)
    ϵ = 1e-4
    nbins = 500
    disc_axis = linspace(0.,quantile(Gamma(T,1.),1.-ϵ),nbins)
    gammas = pdf(Gamma(T,1), disc_axis)
    # compute first two derivatives of 1-cdf as fct of threshold
    dsurvdT = surv(dlgammadT.(T,disc_axis).*gammas,disc_axis)
    dsurvdT2 = surv(dlgammadT2.(T,disc_axis).*gammas,disc_axis)
    # initialize gradient and hessian
    grad = zeros(1,1+size(X,2))
    hes = zeros(1+size(X,2),1+size(X,2))
    for i in 1:length(eff_t0)
        addgrad!(grad, T, eff_t0[i],cs[i], X[i:i,:],dsurvdT,disc_axis)
        addhes!(hes, T, eff_t0[i], cs[i], X[i:i,:],dsurvdT,dsurvdT2,disc_axis)
    end
    return -ll, -grad, -hes
end

#= Given by the user!!!
l(p,t)
grad_small(p,t)
hes_small(p,t)

Step 1) Get:
ll(p,t)
ll(p,t1,t2)
grad_small(p,t1,t2)
hes_small(p,t1,t2)


[T,t] = g(T, β...)

@
grad_big = grad()

#= First STEP:
Go from dP(T,t) to
=#

# Second Step: get hes_big, grad_big

grad_big(T, t1, t2) = hes_()






#=
dlgammadT(T,x) = (x<= 0. ? 0. : log(x))-polygamma(0,T)
dlgammadT2(T,x) = -polygamma(1,T)

grad_lgamma(T,x) = [dlgammadT(T,x) (T-1.)/x-1.]
function grad_lgamma(T,x, c,dsurvdT,disc_axis)
    if !c
        return grad_lgamma(T,x)
    else
        return [interpolate(x, dsurvdT,disc_axis) -pdf(Gamma(T,1.),x)]/(1-cdf(Gamma(T,1.),x))
    end
end

hes_lgamma(T,x) = [dlgammadT2(T,x) 1/x; 1/x -(T-1.)/x^2]
function hes_lgamma(T,x, c,dsurvdT,dsurvdT2,disc_axis)
    if !c
        return hes_lgamma(T,x)
    else
        h11 = interpolate(x, dsurvdT2, disc_axis)
        h12,h22 = -pdf(Gamma(T,1.),x)*grad_lgamma(T,x)
        return [h11 h12; h12 h22]/(1.-cdf(Gamma(T,1.),x))
    end
end

function addgrad!(grad, T,x,c,X_i,dsurvdT,disc_axis)
    grad_small = grad_lgamma(T,x,c,dsurvdT,disc_axis)
    grad[1:1,1] += grad_small[1]
    grad[1:1,2:end] += grad_small[2]*x*X_i
    return
end
function addhes!(hes, T,x,c,X_i,dsurvdT,dsurvdT2,disc_axis)
    grad_small = grad_lgamma(T,x,c,dsurvdT,disc_axis)
    hes_small = hes_lgamma(T,x,c,dsurvdT,dsurvdT2,disc_axis)
    hes[1,1] += hes_small[1,1]
    hes[2:end, 2:end] += x*X_i'*X_i*grad_small[2] + x^2*X_i'*X_i*hes_small[2,2]
    hes[2:end,1:1] += x*X_i'*hes_small[1, 2]
    hes[1:1, 2:end] += x*X_i*hes_small[1, 2]
    return
end

#@step get_likelihood(rand(100),rand([false, true],100), rand(100,5), 3, rand(5))
