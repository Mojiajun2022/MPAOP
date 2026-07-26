# --- Example Usage (Illustrative - to be run in user's script) ---
using MPAOP
function SOfunc(x::Vector{Float64})
    return sum(x .^ 2) # Sphere function
end

function MOfunc(x::Vector{Float64})
    # Example: ZDT1 (assuming x elements are between 0 and 1)
    f1 = x[1]
    g = 1.0 + (9.0 / (length(x) - 1)) * sum(x[2:end])
    h = 1.0 - sqrt(f1 / g)
    f2 = g * h
    return [f1, f2]
end

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

arc_p_s, arc_o_s, conv_s = MOMPA(
    fobj=MOfunc, lb=lb_mo, ub=ub_mo,
    SearchAgents_no=SearchAgents_no, Max_iter=Max_iter, num_objectives=num_obj,
    p0_optional=p_empty, parallelism=:mpi,write_csv_log = true,
    disp=true, saveHDF=false, disp_param=true, hdf_filepath="mo_serial.h5", history_save_interval=Max_iter
)
