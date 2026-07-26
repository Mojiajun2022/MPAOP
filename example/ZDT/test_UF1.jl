using MPAOP
using Plots

include("UF1.jl")

println("\n--- Running MOMPA for UF1 ---")

SearchAgents_no = 200
Max_iter = 500 
dim_mo = 20        
num_obj = 2

lb_mo = [0.0; fill(-1.0, dim_mo - 1)]
ub_mo = [1.0; fill(1.0, dim_mo - 1)]
p_empty = []

archive_prey_uf1, archive_obj_uf1, convergence_uf1 = MPAOP.MOMPA(
    fobj=MOfunc_UF1, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, 
    disp=true,           
    variant = :nmpa,
    disp_param=true,  
    first_stage_ratio=1 / 3,
    second_stage_ratio=2 / 3, 
    P0 = 0.9,
    FADs0 = 0.1,  
    saveHDF=false
)

if !isnothing(archive_obj_uf1) && !isempty(archive_obj_uf1)
    println("\n📊 MOMPA on UF1 found $(size(archive_obj_uf1, 1)) non-dominated solutions. Plotting...")

    f1_obtained = archive_obj_uf1[:, 1]
    f2_obtained = archive_obj_uf1[:, 2]

    f1_true = range(0, 1, length=200)
    f2_true = 1.0 .- sqrt.(f1_true)

    plot(
        f1_true, f2_true, 
        seriestype=:line,
        linewidth=2, 
        color=:red,
        label="True Pareto Front (UF1)",
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on UF1",
        legend=:topright
    )
    
    scatter!(
        f1_obtained, f2_obtained, 
        markersize=3, 
        markerstrokewidth=0,
        color=:blue,
        alpha=0.7,
        label="Obtained Non-Dominated Solutions"
    )
    
    xlims!(0, 1.1)
    ylims!(0, 1.1)

    display(current())
    savefig("mompa_uf1_front.png")
    println("Plot saved as mompa_uf1_front.png")
else
    println("⚠️ MOMPA on UF1 found no non-dominated solutions or archive_obj is nothing.")
end