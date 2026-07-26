# --- Example Usage (Illustrative - to be run in user's script) ---
using MPAOP
function SOfunc(x::Vector{Float64})
    return sum(x.^2) # Sphere function
end

function MOfunc(x::Vector{Float64})
    # Example: ZDT1 (assuming x elements are between 0 and 1)
    f1 = x[1]
    g = 1.0 + (9.0 / (length(x) - 1)) * sum(x[2:end])
    h = 1.0 - sqrt(f1 / g)
    f2 = g * h
    return [f1, f2]
end

function run_examples()
    SearchAgents_no = 50
    Max_iter = 100
    dim_so = 10
    dim_mo = 30
    lb_so = fill(-10.0, dim_so)
    ub_so = fill(10.0, dim_so)
    lb_mo = zeros(dim_mo)
    ub_mo = ones(dim_mo)
    num_obj = 2

    p_empty = []

    println("\n--- Running Single-Objective MPA (Serial, Standard) ---")
    fit_s, pos_s, curve_s = SOMPA(
        fobj=SOfunc, lb=lb_so, ub=ub_so,
        SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
        p0_optional=p_empty, variant=:standard_mpa, parallelism=:serial,
        disp=true, saveHDF=true, hdf_filepath="so_serial_std.h5", history_save_interval=Max_iter
    )
    println("Best Fitness (Serial, Standard): $fit_s")

    println("\n--- Running Single-Objective MPA (Threads, NMPA) ---")
    # Ensure JULIA_NUM_THREADS is set, e.g., `export JULIA_NUM_THREADS=4` before starting Julia
    fit_t, pos_t, curve_t = SOMPA(
        fobj=SOfunc, lb=lb_so, ub=ub_so,
        SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
        p0_optional=p_empty, variant=:nmpa, parallelism=:threads,
        disp=true, saveHDF=true, hdf_filepath="so_threads_nmpa.h5", history_save_interval=Max_iter
    )
    println("Best Fitness (Threads, NMPA): $fit_t")

    # For MPI example, this would be run with `mpiexec -n <num_procs> julia your_script.jl`
    # Inside your_script.jl:
    # if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    #     println("\n--- Running Single-Objective MPA (MPI, Standard) ---")
    # end
    # fit_m, pos_m, curve_m = Run_SingleObjectiveMPA(
    #     fobj=my_single_obj_func, lb=lb_so, ub=ub_so,
    #     SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
    #     p0_optional=p_empty, variant=:standard_mpa, parallelism=:mpi,
    #     disp=true, saveHDF=true, hdf_filepath="so_mpi_std.h5", history_save_interval=Max_iter
    # )
    # if MPI.Comm_rank(MPI.COMM_WORLD) == 0
    #     println("Best Fitness (MPI, Standard): $fit_m")
    # end


    println("\n--- Running Multi-Objective MPA (Serial) ---")
    arc_p_s, arc_o_s, conv_s = MOMPA(
        fobj=MOfunc, lb=lb_mo, ub=ub_mo,
        SearchAgents_no=SearchAgents_no, Max_iter=Max_iter, num_objectives=num_obj,
        p0_optional=p_empty, parallelism=:serial,
        disp=true, saveHDF=true, disp_param = true, hdf_filepath="mo_serial.h5", history_save_interval=Max_iter
    )
    if !isnothing(arc_p_s)
        println("Archive Size (Serial): $(size(arc_p_s, 1))")
    end

    # MPI.Finalize() # If MPI was initialized by one of the calls and script is ending
end
run_examples()
# To run the MPI part of the example, you'd structure your main script
# to call Run_SingleObjectiveMPA or Run_MultiObjectiveMPA with parallelism=:mpi
# and then handle MPI.Finalize() at the end of the script.
# For example:
#
# if abspath(PROGRAM_FILE) == @__FILE__
#   # This block executes when the script is run directly
#
#   # Check if running with MPI and initialize if needed by a function
#   # For functions that handle MPI.Init() internally, this might not be needed here.
#
#   # Example for an MPI run (assuming the function is called by all processes)
#   use_mpi = true # Or determine from ARGS
#
#   if use_mpi
#       # Ensure MPI is initialized if functions expect it to be
#       # MPI.Initialized() || MPI.Init() # Functions handle this
#
#       # Example call to an MPI enabled function
#       fit_m, pos_m, curve_m = Run_SingleObjectiveMPA(
#           fobj=my_single_obj_func, lb=lb_so, ub=ub_so, SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
#           parallelism=:mpi # This tells the function to use MPI logic
#       )
#       if MPI.Comm_rank(MPI.COMM_WORLD) == 0
#           println("MPI SO MPA Best Fitness: ", fit_m)
#       end
#
#       MPI.Finalize() # Finalize MPI at the end of the script if it was initialized
#   else
#       run_examples() # Run serial/threaded examples
#   end
# end

