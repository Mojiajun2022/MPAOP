function MOfunc_ZDT3(x::Vector{Float64})
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
    
    val_for_sqrt = f1 / g
    if val_for_sqrt < 0.0 
        h = 1.0 
    else
        h = 1.0 - sqrt(val_for_sqrt) - (val_for_sqrt) * sin(10.0 * pi * f1)
    end
    f2 = g * h
    
    return [f1, f2]
end