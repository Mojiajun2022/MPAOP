

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
using MPI
using Dates

# The Marine Predator Algorithm(MPA) is first proposed by Afshin Faramarzi et al.(Faramarzi, Afshin, et al. Expert systems with applications 152 (2020): 113377.). 
# The original version is the MATLAB version, I adapted the MATLAB version to the Julia version.
# There are some corresponding improvements have been made, such as the introduction of initial values ​​and the ability to call MPI parallel calculations.



function MPA_MPI(SearchAgents_no::Int64, Max_iter::Int64, p0, lb::Vector{Float64}, ub::Vector{Float64}, dim::Int64, fobj::Function; disp::Bool=true, Fixbox::Bool=true, Write::Bool=false, FADs0=0.2, P0=0.5)
    FADs = FADs0
    P = P0
    Top_predator_pos = fill(0.0, 1, dim)
    Top_predator_fit = Inf
    Convergence_curve = fill(0.0, 1, Max_iter)
    stepsize = fill(0.0, SearchAgents_no, dim)
    fitness = fill(Inf, SearchAgents_no)
    Prey = zeros(dim, SearchAgents_no)
    MPI.Initialized() || MPI.Init()
    root = 0
    comm = MPI.COMM_WORLD
    Size = MPI.Comm_size(MPI.COMM_WORLD)
    rank = MPI.Comm_rank(comm)
    na = mod(SearchAgents_no,Size)
    if na == 0

        na = Size
    end
    SearchAgents_no = SearchAgents_no +(Size-na)
    chunk = ceil(Int, SearchAgents_no / Size)

    if rank == root
        MB = zeros(SearchAgents_no, 1)
    end
    if Write &rank == root
            now_time = Dates.now()
            tes = [(text=string(now_time,"  start running!\n","______________________________"),)]
            CSV.write("fitting_process",tes,  append=true, header=false)
    end
    if rank == root
        Xmin = repeat(lb', SearchAgents_no)
        Xmax = repeat(ub', SearchAgents_no)
        time = 0
        Prey_old = copy(Prey)
        fit_old = copy(fitness)
        if isempty(p0)
            Prey = initialization(SearchAgents_no, dim, ub, lb)'
        else
            rnd = (2 * rand(SearchAgents_no, length(ub))) .- 1
            idx = rnd .> 0
            Prey = repeat(p0', SearchAgents_no) .+ rnd .* (ub .- p0)' .* idx .+ rnd .* (p0 - lb)' .* .!idx
            Prey[1, :] = p0
        end
    end
    Iter = 0
    Prey = MPI.bcast(Prey, root, comm)
    if rank == root
    sign = 1
    end
    while Iter < Max_iter
        GC.gc()
        time0 = @elapsed begin
            
            if rank == root
                if Fixbox
                    Prey = clamp.(Prey, lb', ub')
                end
            end
            Prey = MPI.bcast(Prey, root, comm)
            for i = 1:chunk
                send_buf = nothing
                if rank == 0
                    send_buf = collect((i-1)*Size+1:((i-1)*Size+Size))
                end
                v = MPI.Scatter(send_buf, Int, comm; root=0)

                fa = fobj(Prey[v, :])
                recv_buf = MPI.Gather(fa, comm; root=0)
                if rank == root
                    MMA = copy(recv_buf)
                    MB[(i-1)*Size+1:(i-1)*Size+length(send_buf)] .= vec(MMA)
                end
            end
            if rank == root

                fitness[:] = copy(MB[:])
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
                # println(Indx)
                Prey = Indx' .* Prey_old .+ (iszero.(Indx')) .* Prey
                fitness = Inx .* fit_old .+ (iszero.(Inx)) .* fitness
                fit_old = copy(fitness)
                Prey_old = copy(Prey)
                Elite = repeat(Top_predator_pos', SearchAgents_no)
                CF = (1 - Iter / Max_iter)^(2 * Iter / Max_iter)
                RL = 0.05 * levy(SearchAgents_no, dim, 1.5)
                RB = randn(SearchAgents_no, dim)
                # println(size(Prey))
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
            end
            Prey = MPI.bcast(Prey, root, comm)
            for i = 1:chunk
                send_buf = nothing
                if rank == 0
                    send_buf = collect((i-1)*Size+1:((i-1)*Size+Size))
                end
                v = MPI.Scatter(send_buf, Int, comm; root=0)
                fa = fobj(Prey[v, :])
                recv_buf = MPI.Gather(fa, comm; root=0)
                if rank == root
                    MMA = copy(recv_buf)
                    MB[(i-1)*Size+1:(i-1)*Size+length(send_buf)] .= vec(MMA)
                end
            end
            if rank == root
                fitness[:] = copy(MB[:])
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
                    Prey += CF * (Xmin + rand(SearchAgents_no, dim) .* (Xmax - Xmin)) .* U
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
                # 将Prey广播给所有核

            end
            Prey = MPI.bcast(Prey, root, comm)
            Iter = MPI.bcast(Iter, root, comm)
        end
        if rank == 0
            if disp
                println(string("The time taken for this iteration is", time0, '.', "The total running time is", time, '.', "The estimated remaining time is:", time ./ Iter * (Max_iter - Iter)))
                println("----------------------------------------------")
            end
            time += time0
            if Write
                # CSV.write("fitting_process_chi2", DataFrame(Fit_list))
                # CSV.write("fitting_process_par_result",DataFrame(POS_list))
                tes = [(text=string("As of the", Iter, "-th iteration, the best error is", round(copy(Top_predator_fit), digits=3), "The best parameters are", round.(copy(vec(Top_predator_pos)), digits=3)),)]
                CSV.write("fitting_process", tes, append=true, header=false)
                # text = [
                # (Name = fitting_parm = POS_list), 
                # (fitting_chi2 = Fit_list)
                # ] 
                # CSV.write("fitting_process", text)    
                # writedlm("fitting_process_chi2", Fit_list,delim='\t')
                # writedlm("fitting_process_par_result", POS_list,delim='\t')
            end
        end

        # MPI.Barrier(comm)
    end
    if rank == root
        # round(0.123456, digits = 3)
        Top_predator_pos = round.(copy(vec(Top_predator_pos)), digits=3)
        Top_predator_fit = round.(copy(Top_predator_fit), digits=3)
    end
    
    if Write &rank == root
            end_time = Dates.now()
            tes = [(text=string(end_time,"  finish!\n","______________________________"),)]
            CSV.write("fitting_process",tes, append=true, header=false)
    end
    return Top_predator_fit, Top_predator_pos, Convergence_curve
end