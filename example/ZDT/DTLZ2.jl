function MOfunc_DTLZ2(x::Vector{Float64}, M::Int = 3) # M is the number of objectives
    n = length(x)
    k = n - M + 1 # Number of variables for g function
    
    if n < M
        # Not enough variables for the specified number of objectives
        return fill(Inf, M)
    end

    objectives = zeros(Float64, M)
    
    g = 0.0
    if k > 0 # Check if there are variables for g
        # x_M are the last k variables: x[M], x[M+1], ..., x[n]
        for i = M:n # Indices in x for g function variables
            g += (x[i] - 0.5)^2
        end
    end
    
    common_term = 1.0 + g
    
    for i = 1:M
        product_cos = 1.0
        if i > 1
            for j = 1:(i - 1)
                product_cos *= cos(x[j] * 0.5 * pi)
            end
        end
        
        if i < M
            objectives[i] = common_term * product_cos * sin(x[i] * 0.5 * pi)
        else # Last objective f_M
            objectives[M] = common_term * product_cos 
        end
    end

     for i = 1:M
        term = common_term
        for j = 1:(M - i)
            term *= cos(x[j] * pi / 2.0)
        end
        if i > 1 
            term *= sin(x[M - i + 1] * pi / 2.0)
        end
    end

    fill!(objectives, common_term) 
    for i = 1:M
        for j = 1:(M - i) 
            objectives[i] *= cos(x[j] * 0.5 * pi)
        end
        if i > 1 
        end
    end
    
    
    g_val = 0.0
    
    if n >= M 
        for i_g = M:n
            g_val += (x[i_g] - 0.5)^2
        end
    end
    
    one_plus_g = 1.0 + g_val

    for i_obj = 1:M
        if i_obj == M 
            objectives[M] = one_plus_g * sin(x[1] * 0.5 * pi)
        else
            term = one_plus_g
            for j_cos = 1:(M - i_obj)
                term *= cos(x[j_cos] * 0.5 * pi)
            end
            term *= sin(x[M - i_obj + 1] * 0.5 * pi)
            objectives[i_obj] = term
        end
    end
      fill!(objectives, one_plus_g) 

    for i_obj = 1:M 
        for j_var = 1:(M - i_obj) 
            objectives[i_obj] *= cos(x[j_var] * 0.5 * pi)
        end
        if i_obj < M 
            objectives[i_obj] *= sin(x[M - i_obj] * 0.5 * pi)
        elseif i_obj == M && M > 1 
        end
    end
    
    
    objectives[1] = one_plus_g
    for i = 1:(M-1) 
        objectives[1] *= cos(x[i] * 0.5 * pi)
    end

    for m = 2:M 
        objectives[m] = one_plus_g
        for i = 1:(M-m) 
            objectives[m] *= cos(x[i] * 0.5 * pi)
        end
        objectives[m] *= sin(x[M-m+1] * 0.5 * pi)
    end


    return objectives
end