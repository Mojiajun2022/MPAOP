
using MPAOP 
using Plots
include("ZDT3.jl")


SearchAgents_no = 100 
Max_iter = 250      
dim_mo = 50         
num_obj = 2

lb_mo = zeros(dim_mo) 
ub_mo = ones(dim_mo)  

p_empty = [] 


archive_prey_zdt3, archive_obj_zdt3, convergence_zdt3 = MPAOP.MOMPA(
    fobj=MOfunc_ZDT3, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, 
    disp=true,           
    disp_param=false,     
    saveHDF=true, 
    P0 = 0.5,      
    FADs0 = 0.2,
    hdf_filepath="mo_serial_zdt3.h5", 
    history_save_interval=Max_iter,
    write_csv_log = false
)

if !isnothing(archive_obj_zdt3) && !isempty(archive_obj_zdt3)
    println("\n📊 MOMPA在ZDT3上找到了 $(size(archive_obj_zdt3, 1)) 个非支配解。正在绘图...")

    f1_obtained = archive_obj_zdt3[:, 1]
    f2_obtained = archive_obj_zdt3[:, 2]

    f1_true_values = range(0, 1, length=1000)
    f2_true_values = map(f1_val -> 1.0 - sqrt(f1_val) - f1_val * sin(10.0 * pi * f1_val), f1_true_values)
    plot(
        f1_true_values, f2_true_values, 
        seriestype=:line,
        linewidth=2, 
        color=:red,
        label="True Pareto Front (ZDT3)",
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on ZDT3",
        legend=:bottomleft
    )
    
    scatter!(
        f1_obtained, f2_obtained, 
        markersize=3, 
        markerstrokewidth=0,
        color=:blue,
        alpha=0.7,
        label="Obtained Non-Dominated Solutions"
    )

    display(current())

else
    println("⚠️ MOMPA未能在ZDT3上找到非支配解 (archive_obj_zdt3 为空或nothing)。")
end
