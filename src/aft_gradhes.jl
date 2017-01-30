# # Grad: 1:pdist.M is w.m. of gradlog, pdist.M is -pdf/(1-cdf)
# # Hes: -grad*grad' +
# #  1:pdist.M:1:pdist.M -> (pdist.gradlog(ϕ,aux(x))[i]*pdist.gradlog(ϕ,aux(x))[j]+pdist.heslog(ϕ,aux(x))[i,j]
# # last row: + grad[pdist.M+1]*pdist.gradlog(ϕ,t)
#

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
