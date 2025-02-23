using SymEngine 

@vars x; 

f = x^2 + x/2 - sin(x)/x

diff(f,x) # symbolic diff 

# numerical scheme 
diff_forward(f, x; h = sqrt(eps(Float64))) = (f(x+h) - f(x))/h
diff_central(f, x; h = cbrt(eps(Float64))) = (f(x+h/2) - f(x-h/2))/h
diff_backward(f, x; h = sqrt(eps(Float64))) = (f(x) - f(x-h))/h 

# the esp() method is used to represent a step size between 1.0 and the nearest round number 
diff_complex(f, x; h =1e-20) = imag(f(x+ h*im)) / h 
# imaginary step method 


f = x -> sin(x^2)
x_val = 1.0 
@show diff_forward(f, x_val)
@show diff_central(f, x_val)
@show diff_backward(f, x_val)
@show diff_complex(f, x_val)


struct Dual
    v
    ∂
end 

Base.:+(a::Dual, b::Dual) = Dual(a.v + b.v, a.∂ + b.∂)
Base.:*(a::Dual, b::Dual) = Dual(a.v * b.v, a.v*b.∂ + b.v*a.∂)
Base.log(a::Dual) = Dual(log(a.v), a.∂/a.v)

function Base.max(a::Dual, b::Dual)
    v = max(a.v, b.v)
    ∂ = (a.v > b.v ? a.∂ : a.v < b.v ? b.∂ : NaN)
    return Dual(v, ∂)
end 

function Base.max(a::Dual, b::Int)
    v = max(a.v, b)
    ∂ = a.v > b ? a.∂ : a.v < b ? 0 : NaN
    return Dual(v, ∂)
end 


# bracketing for a unimodal function 
function bracket_minimum(f, x; s = 1e-2, k = 2.0) # s and k are hyperparameters where s is the size of the step and then k is scaling the step
    a, ya = x, f(x)
    b, yb = a + s, f(a + s)
    if yb > ya
        a, b = b, a 
        ya, yb = yb, ya 
        s = -s
    end 
    while true
        c, yc = b + s, f(b + s)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end 
        a, ya, b, yb = b, yb, c, yc
        s *= k
    end 
end 

#functions with no local minima would fail at the above bracket_minimum function 










