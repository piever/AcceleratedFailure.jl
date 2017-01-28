# # Grad: 1:pdist.M is w.m. of gradlog, pdist.M is -pdf/(1-cdf)
# # Hes: -grad*grad' +
# #  1:pdist.M:1:pdist.M -> (pdist.gradlog(ϕ,aux(x))[i]*pdist.gradlog(ϕ,aux(x))[j]+pdist.heslog(ϕ,aux(x))[i,j]
# # last row: + grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#
# aux(x) = gama.quantile(ϕ,(sin(x)+1)/2)
# aux_inverse(t) = asin(2*gama.cdf(ϕ,t)-1)
#
# function cdfgradhessiannew(pdist,ϕ,times)
#     f1 = Fun(x-> 1/2*cos(x)*pdist.gradlog(ϕ,aux(x))[1], -π/2..π/2,50)
#     g1 = cumsum(f1)
#     grad = zeros(pdist.M+1)
#     for t in times
#         grad[1] = (extrapolate(g1,aux_inverse(Inf))-extrapolate(g1,aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#     end
#     grad[2] = -pdist.pdf(ϕ,t)/(1-pdist.cdf(ϕ,t))
#     return grad, f1
# end
#
#
# aux(x) = 1. /(1. -x) - 1.
# aux_(x) = 1. /((1. -x)*(1.-x))
# aux_inverse(t) = 1. - 1. /(1. +t)
#
# function cdfgradhessianvec(pdist,ϕ,ts)
#     f1 = [Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*pdist.gradlog(ϕ,aux(x))[i], aux_inverse(0)..aux_inverse(Inf))
#     for i in 1:pdist.M]
#     g1 = cumsum.(f1)
#     f2 = [Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*
#     (pdist.gradlog(ϕ,aux(x))[i]*pdist.gradlog(ϕ,aux(x))[j]+pdist.heslog(ϕ,aux(x))[i,j]),
#     aux_inverse(0)..aux_inverse(Inf)) for i in 1:pdist.M, j in 1:pdist.M]
#     g2 = cumsum.(f2)
#     grad = zeros(pdist.M+1)
#     hes = zeros(pdist.M+1,pdist.M+1)
#     for t in ts
#         for i in 1:pdist.M
#             grad[i] = (g1[i](aux_inverse(Inf))-g1[i](aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#         end
#         grad[pdist.M+1] = -pdist.pdf(ϕ,t)/(1-pdist.cdf(ϕ,t))
#         hes = -grad*grad'
#         for i in 1:pdist.M
#             for j in 1:pdist.M
#                 hes[i,j] += (g2[i,j](aux_inverse(Inf))-g2[i,j](aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#             end
#         end
#         hes[end, :] += grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#         hes[:,end] = hes[end, :]
#     end
#     return grad, hes
# end
#
#
#
#
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











#
#
# function cdfgradhessiano(pdist,ϕ,t)
#     grad = zeros(length(pdist.M)+1)
#     f1 = Fun(x->pdist.pdf(ϕ,x)*pdist.gradlog(ϕ,x)[1:pdist.M], 0..Inf)
#     g1 = cumsum(f1)
#     grad[1:pdist.M] = (g1(Inf)-g1(t))/(1-pdist.cdf(ϕ,t))
#     grad[pdist.M+1] = -pdist.pdf(ϕ,t)/(1-pdist.cdf(ϕ,t))
#     f2 = Fun(x->pdist.pdf(ϕ,x)*
#     (pdist.gradlog(ϕ,x)[1:pdist.M]*gamma.gradlog(ϕ,x)[1:pdist.M]'+pdist.heslog(ϕ,x)[1:pdist.M,1:pdist.M]),
#     0..Inf)
#     g2 = cumsum(f2)
#     hes = -grad*grad'
#     hes[1:pdist.M, 1:pdist.M] += (g2(Inf)-g2(t))/(1-pdist.cdf(ϕ,t))
#     hes[end, :] += grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#     hes[:,end] = hes[end, :]
#     return grad, hes
# end
#
#
# function cdfgradhes(pdist,ϕ,t)
#     grad = zeros(length(pdist.M)+1)
#     f1 = Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*pdist.gradlog(ϕ,aux(x))[1:pdist.M], aux_inverse(0)..aux_inverse(Inf))
#     g1 = cumsum(f1)
#     grad[1:pdist.M] = (g1(aux_inverse(Inf))-g1(aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#     grad[pdist.M+1] = -pdist.pdf(ϕ,t)/(1-pdist.cdf(ϕ,t))
#     f2 = Fun(x->aux_(x)*pdist.pdf(ϕ,aux(x))*
#     (pdist.gradlog(ϕ,aux(x))[1:pdist.M]*gamma.gradlog(ϕ,aux(x))[1:pdist.M]'+pdist.heslog(ϕ,aux(x))[1:pdist.M,1:pdist.M]),
#     aux_inverse(0)..aux_inverse(Inf))
#     g2 = cumsum(f2)
#     hes = -grad*grad'
#     hes[1:pdist.M, 1:pdist.M] += (g2(aux_inverse(Inf))-g2(aux_inverse(t)))/(1-pdist.cdf(ϕ,t))
#     hes[end, :] += grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#     hes[:,end] = hes[end, :]
#     return grad, hes
# end
