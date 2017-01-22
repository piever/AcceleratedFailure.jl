function newton_raphson(f_,h_!,x; ρ = 0.5, c = 1e-4, tol = 1e-4)
    grad= fill(Inf,length(x))
    hes= ones(length(x),length(x))
    y = 1.
    while norm(grad) > tol
        y = h_!(x, grad, hes)
        search_dir = -(cholfact(Positive, hes)\grad)
        step_size = 1.
        while f_(x+search_dir*step_size) > y+c*step_size*dot(grad',search_dir)
            step_size *= ρ
            step_size > 1e-10 || error("Linesearch failed! Problematic Hessian or Gradient?")
        end
        x = x + search_dir*step_size
    end
    return x, y, grad, hes
end

function newton_raphson(f_, x; ρ = 0.5, c = 1e-4, tol = 1e-4)
    h_! = function(x,grad,hes)
        ForwardDiff.gradient!(grad,f_,x)
        ForwardDiff.hessian!(hes,f_,x)
        return f_(x)
    end
    return newton_raphson(f_, h_!, x; ρ = ρ, c = c, tol = tol)
end
