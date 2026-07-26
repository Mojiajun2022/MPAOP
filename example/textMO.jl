using MPAOP
function my_multi_objective_function(x::Vector{Float64})
    # x is a vector of decision variables
    num_vars = length(x)
    f1 = x[1]

    g = 1.0 + (9.0 / (num_vars - 1)) * sum(x[2:end])
    h = 1.0 - sqrt(f1 / g)
    f2 = g * h

    return [f1, f2] # Must return a vector of objectives
end

SearchAgents_no = 100
Max_iter = 200
dim = 30 # Number of decision variables
lb = zeros(dim) # Lower bounds vector
ub = ones(dim)  # Upper bounds vector
p0 = [] # Optional initial guess for the first agent
num_objectives = 2 # Critical: specify the number of objectives

# To run with MPI:
# mpiexec -n 4 julia your_script_name.jl

pareto_solutions, pareto_objectives, convergence = MOMPA_MPI(
    SearchAgents_no, Max_iter, p0, lb, ub, num_objectives, my_multi_objective_function,
    disp=true, disp_param=true, Fixbox=true, Write=true, saveHDF=true
)
using MPI
if !isnothing(pareto_solutions) && MPI.Comm_rank(MPI.COMM_WORLD) == 0
    println("Optimization finished.")
    println("Number of Pareto optimal solutions found: ", size(pareto_solutions, 1))
    # You can then plot or analyze pareto_solutions and pareto_objectives
    # For example, save them to a CSV:
    # using DataFrames, CSV
    # df_solutions = DataFrame(pareto_solutions, :auto)
    # CSV.write("pareto_solutions.csv", df_solutions)
    # df_objectives = DataFrame(pareto_objectives, :auto)
    # CSV.write("pareto_objectives.csv", df_objectives)
end
