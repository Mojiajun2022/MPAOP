#xlims = [-5,5]
#class dims = 3
function MOfunc_Kursawe(x::Vector{Float64})
    n = length(x)
    if n == 0 || n < 2 
        return [Inf, Inf]
    end

    f1 = 0.0
    for i = 1:(n - 1)
        f1 += -10.0 * exp(-0.2 * sqrt(x[i]^2 + x[i+1]^2))
    end

    f2 = 0.0
    for i = 1:n
        f2 += abs(x[i])^0.8 + 5.0 * sin(x[i]^3)
    end
    
    return [f1, f2]
end