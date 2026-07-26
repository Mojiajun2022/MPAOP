using MPAOP
using Plots # Add this line for plotting
# using Printf # Optional, for more formatted printing if needed

# Your MOfunc (ZDT1) definition - ensure it's correct
function MOfunc(x::Vector{Float64})
    # ZDT1: x_i in [0, 1]
    # Ensure your lb_mo and ub_mo correctly define this range for all dim_mo variables.
    if any(val -> val < 0.0 || val > 1.0, x)
        # Return a very bad objective if x is out of expected [0,1] bounds for ZDT1 variables
        # This helps if the algorithm explores outside, though clamp in MOMPA should handle it.
        # println("Warning: x out of expected [0,1] bounds for ZDT1: ", x) # Optional debug
        return [Inf, Inf] 
    end
    f1 = x[1]
    g = 1.0
    if length(x) > 1
        g += (9.0 / (length(x) - 1)) * sum(x[2:end])
    end
    h = 1.0 - sqrt(f1 / g) # f1/g must be <= 1.0 and >= 0 for sqrt to be real.
                          # g should always be >= 1. f1 is x[1] in [0,1]. So f1/g should be in [0,1].
    f2 = g * h
    return [f1, f2]
end

SearchAgents_no = 100
Max_iter = 100
# dim_so = 10 # Not used for MO example
dim_mo = 100
# lb_so = fill(-10.0, dim_so) # Not used for MO example
# ub_so = fill(10.0, dim_so) # Not used for MO example
lb_mo = zeros(dim_mo) # ZDT1 variables are in [0,1]
ub_mo = ones(dim_mo)  # ZDT1 variables are in [0,1]
num_obj = 2

p_empty = [0,1]

println("🧪 Running MOMPA for ZDT1 visualization...")
# Assuming MOMPA is the correct function name from your MPAOP module
# If you aliased Run_MultiObjectiveMPA to MOMPA, this is fine.
arc_p_s, arc_o_s, conv_s = MOMPA(
    fobj=MOfunc, 
    lb=lb_mo, 
    ub=ub_mo,
    SearchAgents_no=SearchAgents_no, 
    Max_iter=Max_iter, 
    num_objectives=num_obj,
    p0_optional=p_empty, 
    parallelism=:serial, # Keeping serial for straightforward testing
    disp=true,           # Show iteration progress
    disp_param=true,   # Don't need to print params per iteration for this viz
    archive_size_factor=1.0,   
    saveHDF=false,         # Optional: turn off HDF5 saving for quick viz test
    hdf_filepath="mo_serial_zdt1.h5", 
    first_stage_ratio = 1/3,
    second_stage_ratio = 2/3,
    history_save_interval=0
)

if !isnothing(arc_o_s) && !isempty(arc_o_s) # Check if arc_o_s has data
    println("\n📊 Archive found with $(size(arc_o_s, 1)) non-dominated solutions.")
    println("Plotting results...")

    # 1. Extract obtained objective values
    f1_obtained = arc_o_s[:, 1]
    f2_obtained = arc_o_s[:, 2]

    # 2. Generate points for the true Pareto front of ZDT1
    # For ZDT1, the true Pareto front is f2 = 1 - sqrt(f1), where f1 is in [0, 1] and g(x)=1.
    f1_true = range(0, 1, length=200)
    f2_true = 1.0 .- sqrt.(f1_true)

    # 3. Create the plot
    plot(
        f1_true, f2_true, 
        label="True Pareto Front (ZDT1)", 
        linewidth=2, 
        color=:red,
        xlabel="Objective 1 (f1)", 
        ylabel="Objective 2 (f2)",
        title="MOMPA Performance on ZDT1",
        legend=:topright 
    )
    
    scatter!(
        f1_obtained, f2_obtained, 
        label="Obtained Non-Dominated Solutions", 
        markersize=4, 
        markerstrokewidth=0.5,
        color=:blue,
        alpha=0.7 # Transparency
    )

    # Set axis limits for better visualization if needed, e.g., if your results are far off
    # xlims!(0, 1.1) 
    # ylims!(0, 1.1)

    # Display the plot
    # In a script, the plot might show in a new window or in the VS Code plot pane.
    # If running purely from terminal without a display server, you might need to save it.
    # For saving:
    # savefig("mompa_zdt1_front.png") 
    # println("Plot saved as mompa_zdt1_front.png")
    
    # This will attempt to display the plot. 
    # If you are in a REPL or an environment like Jupyter or VS Code with Julia extension, it should display.
    display(current()) # Ensures the plot is shown

else
    println("⚠️ No non-dominated solutions found in the archive (arc_o_s is nothing or empty). Cannot plot.")
end
