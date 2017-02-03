N = 50000
x = randn(N)
y = randn(N)
z = randn(N)
t1 = rand.(Gamma.(10,1*exp.(x-0.3y)))
t2 = rand(Gamma(15,1),N)

W = [(t2[i]>t1[i]) ? Event(t1[i], false) : Event(t2[i], true) for i = 1:N]

df = DataFrame(x = x, y = y, z = z, a = W)

res = aft(a ~ 1 + x +y + x*y+ z, df, PGamma(); tol = 1e-3, c = 1e-4)

true_res = [log(10), 0., 1., -0.3, 0., 0.]
@test_approx_eq_eps res.coefmat.cols[1] true_res 1e-1
println("Aft works:")
println("$(res.coefmat.cols[1])~$(true_res)")
