using MPAOP, Plots
include("DTLZ2.jl")
dim_dtlz2 = 5
num_obj_dtlz2 = 3
lb_dtlz2 = zeros(dim_dtlz2)
ub_dtlz2 = ones(dim_dtlz2)
function generate_dtlz2_true_front_3d(num_points_per_axis::Int = 20)
    x1_vals = range(0, 1, length=num_points_per_axis)
    x2_vals = range(0, 1, length=num_points_per_axis)
    
    f1_true = Float64[]
    f2_true = Float64[]
    f3_true = Float64[]
    
    for x1_val in x1_vals
        for x2_val in x2_vals
            theta1 = x1_val * 0.5 * pi
            theta2 = x2_val * 0.5 * pi

            push!(f1_true, cos(theta1) * cos(theta2))
            push!(f2_true, cos(theta1) * sin(theta2))
            push!(f3_true, sin(theta1))
        end
    end
    return f1_true, f2_true, f3_true
end

println("🧪 Running MOMPA for DTLZ2 (3 Objectives)...")
archive_prey_dtlz2, archive_obj_dtlz2, convergence_dtlz2 = MOMPA(
    fobj=MOfunc_DTLZ2, 
    archive_size_factor = 0.5,
    lb=lb_dtlz2, 
    ub=ub_dtlz2,
    SearchAgents_no=200, 
    Max_iter=500,        
    num_objectives=num_obj_dtlz2,
    parallelism=:serial, 
    variant = :nmpa,
    disp=true
)
  f1_obtained = archive_obj_dtlz2[:, 1]
  f2_obtained = archive_obj_dtlz2[:, 2]
  f3_obtained = archive_obj_dtlz2[:, 3]


  f1_tpf, f2_tpf, f3_tpf = generate_dtlz2_true_front_3d(30) 

#=
  plot3d(
      f1_tpf, f2_tpf, f3_tpf,
      seriestype=:scatter,
      markersize=1.5,
      markerstrokewidth=0,
      markeralpha=0.3, # 透明度
      color=:red,
      label="True Pareto Front (DTLZ2, M=3)",
      xlabel="f1",
      ylabel="f2",
      zlabel="f3",
      title="MOMPA Performance on DTLZ2 (3 Objectives)",
      camera=(30, 30) # 调整视角
  )

  scatter3d!(
      f1_obtained, f2_obtained, f3_obtained,
      markersize=3,
      markerstrokewidth=0.5,
      color=:blue,
      label="Obtained Non-Dominated Solutions",
      alpha=0.8
  )
  
  xlims!(0, 1.1)
  ylims!(0, 1.1)
  zlims!(0, 1.1)

  display(current()) 
  # savefig("mompa_dtlz2_3d_front.png")
  # println("图像已保存为 mompa_dtlz2_3d_front.png")
=#
f1_obtained = archive_obj_dtlz2[:, 1]
f2_obtained = archive_obj_dtlz2[:, 2]
f3_obtained = archive_obj_dtlz2[:, 3]

scatter(
    f1_obtained, 
    f2_obtained, 
    marker_z = f3_obtained, 
    color = :viridis,       
    markersize = 4,
    markerstrokewidth = 0,
    label = "Obtained Solutions (Color = f3)",
    xlabel = "Objective 1 (f1)",
    ylabel = "Objective 2 (f2)",
    title = "MOMPA on DTLZ2 (f1 vs f2, color-coded by f3)",
    colorbar_title = "Objective 3 (f3)", 
    legend = :topright
)


# xlims!(0, 1.1)
# ylims!(0, 1.1)

f1_obtained = archive_obj_dtlz2[:, 1]
f2_obtained = archive_obj_dtlz2[:, 2]


f1_bins = 0:0.05:1.0 
f2_bins = 0:0.05:1.0
using StatsBase
h = fit(Histogram, (f1_obtained, f2_obtained), (f1_bins, f2_bins))

heatmap(
    h.edges[1], 
    h.edges[2], 
    h.weights,  
    seriescolor = :viridis,
    xlabel = "Objective 1 (f1)",
    ylabel = "Objective 2 (f2)",
    title = "Density Heatmap of Solutions in f1-f2 plane (DTLZ2)",
    colorbar_title = "Number of Solutions",
    aspect_ratio = 1
)

display(current())