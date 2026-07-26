# =============================================================================
#  OW.jl -- HDF5 history output / input
# =============================================================================

function _chunk_dims(A::AbstractArray)
    sz = size(A)
    nd = length(sz)
    nd == 0 && return sz
    # chunk along everything but the last (time / solution) axis
    return ntuple(i -> i == nd ? 1 : sz[i], nd)
end

function _write_dset(f, name::String, data, compress::Int)
    if compress > 0 && length(data) > 0
        try
            f[name, chunk=_chunk_dims(data), deflate=compress] = data
            return
        catch
            # fall back to an uncompressed write if the filter is unavailable
        end
    end
    f[name] = data
end

"""
    SaveMPAHistory(filepath, positions_history_slice, performance_history_slice,
                   is_multi_objective;
                   archive_prey_at_save   = nothing,
                   archive_objectives_at_save = nothing,
                   convergence_curve      = nothing,
                   compress               = 0,
                   is_root_process        = true)

Write the population history to `filepath` (HDF5).

Layout (unchanged, so files written by MPAOP ≤ 0.2 stay readable):

    MPA/PositionsHistory   nagents × dim × nsnapshots
    MPA/FitnessHistory     nagents × nsnapshots              (single objective)
    MPA/ObjectivesHistory  nagents × nobj × nsnapshots       (multi objective)
    MPA/ArchivePrey        narchive × dim                    (multi objective)
    MPA/ArchiveObjectives  narchive × nobj                   (multi objective)
    MPA/ConvergenceCurve   niter                             (new, optional)

`compress` (0-9) enables chunked deflate; on long runs it typically shrinks the
file by 3-5× at a negligible CPU cost.  The redundant `rm` + `sleep(0.05)` of
the previous implementation is gone -- `h5open(..., "w")` truncates by itself,
which removes a 50 ms stall from every checkpoint.
"""
function SaveMPAHistory(
    filepath::String,
    positions_history_slice::AbstractArray{Float64,3},
    performance_history_slice,
    is_multi_objective::Bool;
    archive_prey_at_save::Union{AbstractMatrix{Float64},Nothing}=nothing,
    archive_objectives_at_save::Union{AbstractMatrix{Float64},Nothing}=nothing,
    convergence_curve::Union{AbstractVector{Float64},Nothing}=nothing,
    compress::Int=0,
    is_root_process::Bool=true
)
    is_root_process || return nothing

    folder_path = dirname(filepath)
    if !isempty(folder_path) && !isdir(folder_path)
        try
            mkpath(folder_path)
        catch e
            @warn "Could not create directory $folder_path" exception = e
            return nothing
        end
    end

    try
        h5open(filepath, "w") do f
            _write_dset(f, "MPA/PositionsHistory", Array(positions_history_slice), compress)
            if is_multi_objective
                _write_dset(f, "MPA/ObjectivesHistory", Array(performance_history_slice), compress)
                if archive_prey_at_save !== nothing && size(archive_prey_at_save, 1) > 0
                    _write_dset(f, "MPA/ArchivePrey", Array(archive_prey_at_save), compress)
                end
                if archive_objectives_at_save !== nothing && size(archive_objectives_at_save, 1) > 0
                    _write_dset(f, "MPA/ArchiveObjectives", Array(archive_objectives_at_save), compress)
                end
            else
                _write_dset(f, "MPA/FitnessHistory", Array(performance_history_slice), compress)
                if archive_prey_at_save !== nothing && size(archive_prey_at_save, 1) > 0
                    _write_dset(f, "MPA/BestPosition", Array(archive_prey_at_save), compress)
                end
            end
            if convergence_curve !== nothing && length(convergence_curve) > 0
                f["MPA/ConvergenceCurve"] = Array(convergence_curve)
            end
            try
                HDF5.attributes(f)["created"] = string(Dates.now())
                HDF5.attributes(f)["package"] = "MPAOP.jl"
            catch
            end
        end
    catch h5_error
        @error "HDF5 write failed for $filepath" exception = h5_error
        rethrow(h5_error)
    end
    return nothing
end

"""
    ReadMPAHistory(filepath)
        -> positions_history, performance_history, is_multi_objective,
           archive_prey, archive_objectives

Read back a file written by [`SaveMPAHistory`](@ref).  Returns `nothing` for
whatever is absent.
"""
function ReadMPAHistory(filepath::String)
    if !isfile(filepath)
        @warn "File not found: $filepath"
        return nothing, nothing, false, nothing, nothing
    end

    positions_history = nothing
    performance_history = nothing
    is_multi_objective_data = false
    archive_prey = nothing
    archive_objectives = nothing

    try
        h5open(filepath, "r") do file_handle
            if haskey(file_handle, "MPA/PositionsHistory")
                positions_history = read(file_handle["MPA/PositionsHistory"])
            else
                @warn "Dataset 'MPA/PositionsHistory' not found in $filepath"
            end

            if haskey(file_handle, "MPA/ObjectivesHistory")
                performance_history = read(file_handle["MPA/ObjectivesHistory"])
                is_multi_objective_data = true
                haskey(file_handle, "MPA/ArchivePrey") &&
                    (archive_prey = read(file_handle["MPA/ArchivePrey"]))
                haskey(file_handle, "MPA/ArchiveObjectives") &&
                    (archive_objectives = read(file_handle["MPA/ArchiveObjectives"]))
            elseif haskey(file_handle, "MPA/FitnessHistory")
                performance_history = read(file_handle["MPA/FitnessHistory"])
                is_multi_objective_data = false
            else
                @warn "Neither 'MPA/ObjectivesHistory' nor 'MPA/FitnessHistory' found in $filepath"
            end
        end
    catch e
        @error "Error reading HDF5 file $filepath" exception = e
        return nothing, nothing, false, nothing, nothing
    end

    return positions_history, performance_history, is_multi_objective_data,
    archive_prey, archive_objectives
end

"""
    ReadMPAConvergence(filepath) -> Union{Vector{Float64}, Nothing}

Convenience reader for the optional `MPA/ConvergenceCurve` dataset.
"""
function ReadMPAConvergence(filepath::String)
    isfile(filepath) || return nothing
    curve = nothing
    h5open(filepath, "r") do f
        haskey(f, "MPA/ConvergenceCurve") && (curve = read(f["MPA/ConvergenceCurve"]))
    end
    return curve
end
