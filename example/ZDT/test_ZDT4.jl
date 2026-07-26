using MPAOP
using Plots
# best for test nmpa
include("ZDT4.jl") 

println("\n--- Running MOMPA for ZDT4 ---")

SearchAgents_no = 200
Max_iter = 100      
dim_mo = 50         
num_obj = 2

lb_mo = [0.0; fill(-5.0, dim_mo - 1)]
ub_mo = [1.0; fill(5.0, dim_mo - 1)]
p_empty = []

archive_prey_zdt4, archive_obj_zdt4, convergence_zdt4 = MPAOP.MOMPA(
    fobj=MOfunc_ZDT4, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, 
    archive_size_factor= 0.5,
    variant = :nmpa,
    disp=true,           
    disp_param=true,     
    saveHDF=false
)

if !isnothing(archive_obj_zdt4) && !isempty(archive_obj_zdt4)
    println("\n📊 MOMPA on ZDT4 found $(size(archive_obj_zdt4, 1)) non-dominated solutions. Plotting...")

    f1_obtained = archive_obj_zdt4[:, 1]
    f2_obtained = archive_obj_zdt4[:, 2]

    f1_true = range(0, 1, length=200)
    f2_true = 1.0 .- sqrt.(f1_true)

    plot(
        f1_true, f2_true, 
        seriestype=:line,
        linewidth=2, 
        color=:red,
        label="True Pareto Front (ZDT4, g(x)=1)",
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on ZDT4",
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
    


    display(current())
    savefig("mompa_zdt4_front.png")
    println("Plot saved as mompa_zdt4_front.png")
else
    println("⚠️ MOMPA on ZDT4 found no non-dominated solutions or archive_obj is nothing.")
end