# --- Example Usage (Illustrative - to be run in user's script) ---
using MPAOP
function fobj(x) # The function to be fitted
    f1 = abs(x[1] + x[2]) - abs(x[3])
    f2 = x[1] * x[2] * x[3] + 18
    f3 = x[1]^2 * x[2] + 3 * x[3]
    f = abs(f1) + abs(f2) + abs(f3)
    return f
end
SearchAgents_no = 200
Max_iter = 100
dim_so = 3
dim_mo = 30
lb_so = fill(-10.0, dim_so)
ub_so = fill(10.0, dim_so)
lb_mo = zeros(dim_mo)
ub_mo = ones(dim_mo)
num_obj = 2

p_empty = []

fit_s, pos_s, curve_s = MOMPA(
    fobj=fobj, lb=lb_so, ub=ub_so,
    SearchAgents_no=SearchAgents_no, Max_iter=Max_iter,
    p0_optional=p_empty, variant=:standard_mpa, parallelism=:serial,use_gaussian_perturbation = false,gaussian_perturb_on_all_dims = false,
    disp=true, disp_param=true, saveHDF=true, hdf_filepath="so_serial_std.h5", history_save_interval=Max_iter
)



