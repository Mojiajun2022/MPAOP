"""
    MPAOP

Marine Predators Algorithm (MPA) for Julia -- single and multi objective,
serial / multithreaded / MPI / hybrid.

Main entry points: [`SOMPA`](@ref) and [`MOMPA`](@ref).

Original algorithm: A. Faramarzi, M. Heidarinejad, S. Mirjalili, A. H. Gandomi,
*Marine Predators Algorithm: A nature-inspired metaheuristic*,
Expert Systems with Applications 152 (2020) 113377.
"""
module MPAOP

using Random
using LinearAlgebra
using SpecialFunctions
using Dates
using CSV
using HDF5
using MPI

# low level kernels and helpers
include("kernels.jl")
export levy, levy!

# population initialisation
include("initset.jl")
export initialization

# evaluation back-ends (serial / threads / MPI / hybrid)
include("evaluate.jl")

# multi-objective machinery
include("MOfun.jl")
export non_dominated_sort, calculate_crowding_distance!, update_archive,
    select_elite_from_archive, is_dominated

# quality indicators
include("metrics.jl")
export hypervolume, igd, gd, spacing_metric, max_spread, pareto_filter

# HDF5 history I/O
include("OW.jl")
export SaveMPAHistory, ReadMPAHistory, ReadMPAConvergence

# drivers
include("MPA.jl")
export SOMPA, MOMPA

# MPAOP 0.1 positional API (thin wrappers over SOMPA)
include("legacy.jl")
export MPA, MPA_MPI, confidence_interval

end # module
