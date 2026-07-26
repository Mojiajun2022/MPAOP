#xlims = [0,1]
#class dims = 10
function MOfunc_ZDT6(x::Vector{Float64})
    n = length(x)
    if n == 0
        return [Inf, Inf]
    end
    if n < 2 
    end

    f1 = 1.0 - exp(-4.0 * x[1]) * (sin(6.0 * pi * x[1]))^6
    
    g = 0.0
    if n > 1
        sum_xi_for_g = 0.0
        for i = 2:n
            sum_xi_for_g += x[i]
        end
        g = 1.0 + 9.0 * (sum_xi_for_g / (n - 1))^0.25
    else # n == 1
        g = 1.0 # g(x)=1 if only x1 exists
    end
    
    val_for_f1_g_sq = f1 / g
    if val_for_f1_g_sq < 0.0 || g <= 0.0 
        f2 = Inf
    else
        h = 1.0 - (val_for_f1_g_sq)^2
        f2 = g * h
    end
        
    return [f1, f2]
end