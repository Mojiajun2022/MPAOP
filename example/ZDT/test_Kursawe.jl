using MPAOP
using Plots
include("Kursawe.jl") 

println("\n--- Running MOMPA for Kursawe function ---")

SearchAgents_no = 100    
Max_iter = 300          
dim_mo = 3              
num_obj = 2             

lb_mo = fill(-5.0, dim_mo)  
ub_mo = fill(5.0, dim_mo)   
p_empty = []

archive_prey_kursawe, archive_obj_kursawe, convergence_kursawe = MPAOP.MOMPA(
    fobj=MOfunc_Kursawe, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, 
    disp=true,           
    disp_param=false,     
    saveHDF=false        
)

if !isnothing(archive_obj_kursawe) && !isempty(archive_obj_kursawe)
    num_solutions_found = size(archive_obj_kursawe, 1)
    println("\n📊 MOMPA on Kursawe found $(num_solutions_found) non-dominated solutions. Plotting...")

    f1_obtained = archive_obj_kursawe[:, 1]
    f2_obtained = archive_obj_kursawe[:, 2]
    
    plot( 
        f1_obtained, f2_obtained, 
        seriestype=:scatter,
        markersize=3, 
        markerstrokewidth=0,
        color=:blue,
        alpha=0.7,
        label="Obtained Non-Dominated Solutions",
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on Kursawe Function (n=$(dim_mo))",
        legend=:best 
    )
    
    println("\n--- How to Interpret Kursawe Function Results ---")
    println("1. The True Pareto Front (TPF) of the Kursawe function is known to be disconnected, consisting of several separate segments.")
    println("2. Visually compare your plot with known TPF shapes for Kursawe (n=$(dim_mo)) from literature or benchmark collections (e.g., from PlatEMO, Pymoo, or MOEA research papers).")
    println("   - Do your solutions appear to form similar disconnected segments?")
    println("   - Are the objective values in the expected ranges? (f1 typically between approx. -20 and -3 for n=3).")
    println("3. For quantitative assessment, you would need a reference set of true Pareto front points for Kursawe to calculate metrics like IGD (Inverted Generational Distance).")

    display(current()) 
else
    println("⚠️ MOMPA on Kursawe found no non-dominated solutions or archive_obj is nothing.")
end
