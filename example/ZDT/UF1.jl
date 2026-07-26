# suggest dim_mo = 30
# lb = [0.0; fill(-1.0, 29)]
# ub = [1.0; fill(1.0, 29)]

function MOfunc_UF1(x::Vector{Float64})
    n = length(x)
    if n == 0
        return [Inf, Inf]
    end

    count1 = 0
    sum1 = 0.0
    count2 = 0
    sum2 = 0.0

    for j = 2:n
        yj = x[j] - sin(6.0 * pi * x[1] + (j * pi / n))
        yj_sq = yj^2
        if j % 2 == 1 
            sum1 += yj_sq
            count1 += 1
        else 
            sum2 += yj_sq
            count2 += 1
        end
    end

    f1 = x[1]
    if count1 > 0
        f1 += (2.0 / count1) * sum1
    end

    f2 = 1.0 - sqrt(x[1]) 
    if count2 > 0
        f2 += (2.0 / count2) * sum2
    end
    
    return [f1, f2]
end