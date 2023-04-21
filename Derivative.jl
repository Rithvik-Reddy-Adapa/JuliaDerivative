"""
Derivative(f, x::Real; h::Real = 1e-8, n::Int = 1)
`n^th` order derivative of a function `f` at `x`. `h` is delta change
"""
function Derivative(f, x::Real; h::Real = 1e-8, n::Int = 1)
    x = BigFloat(x)
    h = BigFloat(h)
    n = Int(round(n))
    if n == 0
        dfdx = f(x)
    elseif n == 1
        A = x + h
        B = x - h
        dfdx = (f(A) - f(B)) / (2*h)
    else
        A = x + h
        B = x - h
        dfdx = (Derivative(f, A, h = h, n = n-1) - Derivative(f, B, h = h, n = n-1)) / (2*h)
    end
    return dfdx
end


"""
NdDerivative(f, x; n, h = nothing, m = nothing)
Partial derivative of function `f` with `m` input variables at `x`. `n` is order of each variable as a vector. `h` is delta change.
"""
function NdDerivative(f, x; n, h = nothing, m = nothing)
    x = BigFloat.(x)
    m = (m == nothing ? length(x) : m)
    h = (h == nothing ? zeros(length(x)).+1e-8 : h)
    h = BigFloat.(h)
    if m == 0
        dfdx = BigFloat.(f(x...))
    elseif n[m] == 0
        dfdx = NdDerivative(f, x, n = copy(n), m = m-1, h = h)
    elseif n[m] == 1
        A = copy(x)
        A[m] = A[m] + h[m]
        B = copy(x)
        B[m] = B[m] - h[m]
        n[m] = 0
        dfdx = (NdDerivative(f, A, n = copy(n), m = m-1, h = h) - NdDerivative(f, B, n = copy(n), m = m-1, h = h)) / (2*h[m])
    else
        A = copy(x)
        A[m] = A[m] + h[m]
        B = copy(x)
        B[m] = B[m] - h[m]
        n[m] = n[m] - 1
        dfdx = (NdDerivative(f, A, n = copy(n), m = m, h = h) - NdDerivative(f, B, n = copy(n), m = m, h = h)) / (2*h[m])
    end
    return dfdx
end

#=
function Derivative1(f, x; h = 1e-5, n = 1)
    if n == 1
        A = x + (h/2)
        B = x - (h/2)
        dfdx = (f(A) - f(B)) / h
    else
        A = x + (h/2^(n-1))
        B = x - (h/2^(n-1))
        dfdx = (Derivative(f, A, h = h/2, n = n-1) - Derivative(f, B, h = h/2, n = n-1)) / h
    end
    return dfdx
end
=#
;
