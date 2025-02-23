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




