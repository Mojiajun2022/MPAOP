# suggest dim_mo = 10
#lb = [0.0; fill(-5.0, dim_mo-1)]
#ub = [1.0; fill(5.0, dim_mo-1)]
function MOfunc_ZDT4(x::Vector{Float64})
    n = length(x)
    if n == 0
        return [Inf, Inf]
    end
    if n < 2 
    end

    f1 = x[1]
    
    g = 1.0 + 10.0 * (n - 1)
    sum_terms_g = 0.0
    for i = 2:n
        sum_terms_g += x[i]^2 - 10.0 * cos(4.0 * pi * x[i])
    end
    g += sum_terms_g
    
    h_val_inside_sqrt = f1 / g
    if h_val_inside_sqrt < 0.0 || g <= 0.0 
        f2 = Inf 
    else
        h = 1.0 - sqrt(h_val_inside_sqrt)
        f2 = g * h
    end
    
    return [f1, f2]
end