function accelerated_life(times, censored, X, T, r, β)
    R = X*β
    for i = 1:size(times)
        if censored[i]
            nll -= log(1-cdf(Gamma(T, 1/R[i]), times[i]))
            gradient -=

        else
            nll -= pdf(Gamma(T, 1/R[i]), times[i])
            gradient
end

gradgamma(T,x) = [log(x)-polygamma(0,T) (T-1.)/x-1. ]
hessgamma(T,x) = [-polygamma(1,T) 1/x; 1/x -(T-1.)/x^2]

function gradgamma(T,x, censored)
    if !censored
        return gradgamma(T,x)
    else
        return [integra(x, t-> gradgamma(T,t)[1]) -Gamma(T,x)]/(1-cdf(Gamma(T,1),x))
    end
end

function hessgamma(T,x, censored)
    if !censored
        return gradgamma(T,x)
    else
        h11 = integra(x, t-> hessgamma(T,t)[1,1])
        h12,h22 = -Gamma(T,x)*gradgamma(T,x)
        return [h11 h12; h21 h22]/(1-cdf(Gamma(T,1),x))
    end
end

# gradTxofb(T,β) = [1., x*X_i...]
# function hessTxofb(T,β)
#     hes = zeros(lenght(X_i)+1,lenght(X_i)+1)
#     hes[1,1] = 1.
#     hes[2:end, 2:end] = x*X_i*X_i'
#     return hes
# end

gradbeta(T,x,X_i, β,gradpdf) = hcat(gradpdf(T,x)[1], gradpdf(T,x)[2]*x*X_i)
function hessbeta(T,x,X_i,β,gradpdf,hespdf)
    hes = zeros(lenght(X_i)+1,lenght(X_i)+1)
    hes[1,1] = hespdf[1,1]
    hes[2:end, 2:end] = (x*X_i'*X_i)*gradpdf[2] + x^2*X_i'*X_i*gradpdf[2]
    hes[2:end,1] = x*X_i'*hesspdf[1, 2]
    hes[1, 2:end] = hes[2:end,1]'
    return hes
end

function get_likelihood(times, censored, X, T, β)
    eff_times = times.*exp(-X*β)
    ll = sum(log.(pdf(Gamma(T,1),eff_times[!censored])))+
         sum(log.(1-cdf(Gamma(T,1),eff_times[censored])))
    dpdx =
