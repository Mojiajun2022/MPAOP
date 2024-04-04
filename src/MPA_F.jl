using ThreadPools
using Distributed
using Random
using SpecialFunctions
using Base
using FiniteDiff
using PositiveFactorizations
using LinearAlgebra
using Distributions
using CSV


# The Marine Predator Algorithm(MPA) is first proposed by Afshin Faramarzi et al.(Faramarzi, Afshin, et al. Expert systems with applications 152 (2020): 113377.). 
# The original version is the MATLAB version, I adapted the MATLAB version to the Julia version.
# There are some corresponding improvements have been made, such as the introduction of initial values ​​and the ability to call MPI parallel calculations.


function initialization(SearchAgents_no::Int64, dim::Int64, ub::Vector{Float64}, lb::Vector{Float64})
    Positions = zeros((SearchAgents_no), dim)
    Boundary_no = size(ub, 2)
    ub = ub'
    lb = lb'
    if Boundary_no == 1
        Positions = rand(SearchAgents_no, dim) .* (ub .- lb) .+ lb
    end
    if Boundary_no > 1
        for i = 1:dim
            ub_i = ub[i]
            lb_i = lb[i]
            Positions[:, i] = rand(SearchAgents_no) .* (ub_i - lb_i) .+ lb_i
        end
    end
    Positions = Positions'
    return Positions
end



function levy(n::Int64, m::Int64, beta::Float64)
    num = gamma(1 + beta) * sin(pi * beta / 2)
    den = gamma((1 + beta) / 2) * beta * 2^((beta - 1) / 2)
    sigma_u = (num / den)^(1 / beta)

    u = randn(n, m) .* sigma_u
    v = randn(n, m)
    z = u ./ abs.(v) .^ (1 / beta)
    return z
end


function confidence_interval(param::Vector{Float64}, func::Function, level::Float64)
    # 计算Jacobian矩阵
    hessian_fit = FiniteDiff.finite_difference_hessian(func, param)
    #J = Zygote.jacobian(func, param)
    Vchol = Matrix(cholesky(Positive, inv(hessian_fit)).U)
    # 计算协方差矩阵
    pvar = Vchol' * Vchol
    std_fit = sqrt.(diag(pvar))

    # 假设参数服从正态分布，根据置信水平计算临界值
    z = quantile(Normal(), (1 + level) / 2)
    # 计算置信区间
    ci = zeros(length(param), 2)
    ci[:, 1] = param .+ z * std_fit * -1
    ci[:, 2] = param .+ z * std_fit * 1
    # 返回置信区间
    return ci
end

function MPA(SearchAgents_no::Int64, Max_iter::Int64, p0, lb::Vector{Float64}, ub::Vector{Float64}, dim::Int64, fobj::Function; disp::Bool=true, Fixbox::Bool=true, Threads_parallel::Bool=false, Write::Bool=false, FADs0=0.2, P0=0.5)
    # const FADs = 0.2 #The bigger, the more local
# const P = 0.5 #The smaller, the more local
    FADs = FADs0
    P = P0
    Top_predator_pos = fill(0.0, 1, dim)
    Top_predator_fit = Inf
    Convergence_curve = fill(0.0, 1, Max_iter)
    stepsize = fill(0.0, SearchAgents_no, dim)
    fitness = fill(Inf, SearchAgents_no)
    Prey = zeros(dim, SearchAgents_no)

    if isempty(p0)
        Prey = initialization(SearchAgents_no, dim, ub, lb)'
    else
        rnd = (2 * rand(SearchAgents_no, length(ub))) .- 1
        idx = rnd .> 0
        Prey = repeat(p0', SearchAgents_no) .+ rnd .* (ub .- p0)' .* idx .+ rnd .* (p0 - lb)' .* .!idx
        Prey[1, :] = p0
    end

    Prey_old = copy(Prey)
    Xmin = repeat(lb', SearchAgents_no)
    Xmax = repeat(ub', SearchAgents_no)
    Iter = 0
    fit_old = copy(fitness)
    time = 0

    while Iter < Max_iter
        # GC.gc()
        time0 = @elapsed begin
            if Fixbox
                Prey = clamp.(Prey, lb', ub')
            end
            if Threads_parallel
                ThreadPools.@qthreads for i in 1:size(Prey, 1)
                    fitness[i] = fobj(Prey[i, :])
                end
            else
                for i in 1:size(Prey, 1)
                    fitness[i] = fobj(Prey[i, :])
                end
            end
            min_F, ind1 = findmin(fitness)
            if min_F < Top_predator_fit
                Top_predator_fit = copy(min_F)
                Top_predator_pos = copy(Prey[ind1, :])
            end
            if Iter == 0
                fit_old = copy(fitness)
                Prey_old = copy(Prey)
            end
            Inx = fit_old .< fitness
            Indx = repeat(Inx', dim)
            Prey = Indx' .* Prey_old .+ (iszero.(Indx')) .* Prey
            fitness = Inx .* fit_old .+ (iszero.(Inx)) .* fitness
            fit_old = copy(fitness)
            Prey_old = copy(Prey)
            Elite = repeat(Top_predator_pos', SearchAgents_no,1)
            CF = (1 - Iter / Max_iter)^(2 * Iter / Max_iter)
            RL = 0.05 * levy(SearchAgents_no, dim, 1.5)
            RB = randn(SearchAgents_no, dim)
            # println(Prey)
            for i in 1:size(Prey, 1)
                for j in 1:size(Prey, 2)
                    R::Float64 = rand()
                    if Iter < Max_iter / 3
                        stepsize[i, j]::Float64 = RB[i, j] * (Elite[i, j] - RB[i, j] * Prey[i, j] + rand()-0.5)
                        Prey[i, j] += P * R * stepsize[i, j]
                    elseif Iter > Max_iter / 3 && Iter < 2 * Max_iter / 3
                        if i > size(Prey, 1,) / 2
                            stepsize[i, j]::Float64 = RB[i, j] * (RB[i, j] * Elite[i, j] - Prey[i, j] + rand()-0.5)
                            Prey[i, j] = Elite[i,j]+P * CF * stepsize[i, j]
                        else
                            stepsize[i, j]::Float64 = RL[i, j] * (Elite[i, j] - RL[i, j] * Prey[i, j] + rand()-0.5)
                            Prey[i, j] += P * R * stepsize[i, j]
                        end
                    else
                        stepsize[i, j]::Float64 = RL[i, j] * (RL[i, j] * Elite[i, j] - Prey[i, j] + rand()-0.5)
                        Prey[i, j] = Elite[i,j] + P * CF * stepsize[i, j]
                    end
                end
            end
            if Fixbox
                Prey = clamp.(Prey, lb', ub')
            end
            #Threads.@threads 
           
            if Threads_parallel
                ThreadPools.@qthreads for i in 1:size(Prey, 1)
                    fitness[i] = fobj(Prey[i, :])
                end
            else
                for i in 1:size(Prey, 1)
                    fitness[i] = fobj(Prey[i, :])
                end
            end
            min_F, ind1 = findmin(fitness)
            if min_F < Top_predator_fit
                Top_predator_fit = copy(min_F)
                Top_predator_pos = copy(Prey[ind1, :])
            end
            if Iter == 0
                fit_old = copy(fitness)
                Prey_old = copy(Prey)
            end

            Inx = fit_old .< fitness
            Indx = repeat(Inx', dim)
            Prey = Indx' .* Prey_old .+ (iszero.(Indx')) .* Prey
            fitness = Inx .* fit_old .+ (iszero.(Inx)) .* fitness
            fit_old = copy(fitness)
            Prey_old = copy(Prey)
            
            if rand() < FADs
                U = rand(SearchAgents_no, dim) .< FADs
                Prey += CF * ((Xmin + rand(SearchAgents_no, dim) .* (Xmax - Xmin)) .* U)
            else
                r = rand()
                Rs = size(Prey, 1)
                stepsize = (FADs * (1 - r) + r) * (Prey[randperm(Rs), :] - Prey[randperm(Rs), :])
                Prey += copy(stepsize)
            end
            Iter += 1
            Convergence_curve[Iter] = copy(Top_predator_fit)
            if disp
                println(string("*As of the", Iter, "-th iteration, the best error is", round(copy(Top_predator_fit), digits=3)))
                println(string("The best parameters are", round.(copy(vec(Top_predator_pos)), digits=3)))

            end
        end
        if disp
            println(string("The time taken for this iteration is", time0, '.', "The total running time is", time, '.', "The estimated remaining time is:", time ./ Iter * (Max_iter - Iter)))
            println("----------------------------------------------")
        end
        time += time0
        if Write
            tes = [(text=string("As of the", Iter, "-th iteration, the best error is", round(copy(Top_predator_fit), digits=3), "The best parameters are", round.(copy(vec(Top_predator_pos)), digits=3)),)]
            CSV.write("fitting_process", tes, append=true, header=false)

        end
    end
    Top_predator_pos = round.(copy(vec(Top_predator_pos)), digits=3)
    Top_predator_fit = round.(copy(Top_predator_fit), digits=3)
    level = 0.95
    return Top_predator_fit, Top_predator_pos, Convergence_curve
end

