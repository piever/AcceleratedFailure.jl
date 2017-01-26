F(v) = pdf(Gamma(v[1],1),v[2])
const x = rand(10)
const β = rand(10)/10
M = 1
N = 10
ϕ = 10.
t = 5.
τ = t*exp(-dot(β,x))

G(s) = F([s[1],t*exp(-dot(x,s[2:end]))])

dFdϕτ = Calculus.gradient(F,[ϕ,τ])
d2Fdϕτ2 = Calculus.hessian(F,[ϕ,τ])

grad, hes = Survival.aft_gradhes(dFdϕτ, d2Fdϕτ2,τ,x,M,N,β)

grad1 = Calculus.gradient(s->G(s),vcat(ϕ,β))
hes1 = Calculus.hessian(s->G(s),vcat(ϕ,β))
