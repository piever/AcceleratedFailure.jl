function get_likelihood(S, X, ϕ, β, pdist)
    Xβ = X*β
    eff_S = exp(-Xβ).*S
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
