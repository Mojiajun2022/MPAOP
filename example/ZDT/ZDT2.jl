function MOfunc_ZDT2(x::Vector{Float64})
    n = length(x)
    if n == 0
        return [Inf, Inf]
    end
    f1 = x[1]
    
    g = 0.0
    if n > 1
        g = 1.0 + (9.0 / (n - 1)) * sum(x[2:end])
    else # if n == 1
        g = 1.0 
    end
    
    h = 1.0 - (f1 / g)^2
    f2 = g * h
    
    return [f1, f2]
end