using MPAOP
using Plots

include("ZDT6.jl")

println("\n--- Running MOMPA for ZDT6 ---")

SearchAgents_no = 200
Max_iter = 100 
dim_mo = 20          
num_obj = 2

lb_mo = zeros(dim_mo)
ub_mo = ones(dim_mo)
p_empty = []

archive_prey_zdt6, archive_obj_zdt6, convergence_zdt6 = MPAOP.MOMPA(
    fobj=MOfunc_ZDT6, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, 
    archive_size_factor= 0.5,
    disp=true,   
    variant = :nmpa,        
    disp_param=false,     
    saveHDF=false 
)

if !isnothing(archive_obj_zdt6) && !isempty(archive_obj_zdt6)
    println("\n📊 MOMPA on ZDT6 found $(size(archive_obj_zdt6, 1)) non-dominated solutions. Plotting...")

    f1_obtained = archive_obj_zdt6[:, 1]
    f2_obtained = archive_obj_zdt6[:, 2]
    x1_true_range = range(0, 1, length=1000)
    f1_true = map(x1_val -> 1.0 - exp(-4.0 * x1_val) * (sin(6.0 * pi * x1_val))^6, x1_true_range)
     x1_vals_tpf = range(0, 1, length=500)
    f1_tpf_pts = [1.0 - exp(-4.0 * x1) * (sin(6.0 * pi * x1))^6 for x1 in x1_vals_tpf]
    f2_tpf_pts = 1.0 .- (f1_tpf_pts).^2
    
    sort_indices = sortperm(f1_tpf_pts)
    f1_tpf_sorted = f1_tpf_pts[sort_indices]
    f2_tpf_sorted = f2_tpf_pts[sort_indices]


    plot(
        f1_tpf_sorted, f2_tpf_sorted, 
        seriestype=:line,
        linewidth=2, 
        color=:red,
        label="True Pareto Front (ZDT6, g(x)=1)",
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on ZDT6",
        legend=:bottomright 
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
    savefig("mompa_zdt6_front.png")
    println("Plot saved as mompa_zdt6_front.png")
else
    println("⚠️ MOMPA on ZDT6 found no non-dominated solutions or archive_obj is nothing.")
end