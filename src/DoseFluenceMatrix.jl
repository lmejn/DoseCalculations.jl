export dose_fluence_matrix, dose_fluence_matrix!

#--- Dose-Fluence Matrix-------------------------------------------------------

"""
    dose_fluence_matrix(pos, bixels, bixels, surf, calc)

Compute a dose-fluence matrix from dose positions, bixels, external surface and
dose calculation algorithm.

See `dose_fluence_matrix!` for implementation.
"""
function dose_fluence_matrix(pos, bixels::AbstractVector{<:Bixel}, gantry::GantryPosition,
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)
    D = spzeros(length(pos), length(bixels))
    dose_fluence_matrix!(D, pos, bixels, gantry, surf, calc)
end

function dose_fluence_matrix(pos, beamlets::AbstractVector{<:Beamlet},
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)
    D = spzeros(length(pos), length(beamlets))
    dose_fluence_matrix!(D, pos, beamlets, surf, calc)
end

"""
    dose_fluence_matrix!(D, pos, bixels, bixels, surf, calc)

Compute a dose-fluence matrix from dose positions, bixels, external surface and
dose calculation algorithm.

Requires the `point_kernel!` method to be defined for the given dose calculation
algorithm (`calc`). `point_kernel!` computes the dose calculated from the set of
bixels a given dose point. Stores result in `D`.
"""
function dose_fluence_matrix!(D::SparseMatrixCSC, pos, bixels::AbstractVector{<:Bixel},
                              gantry::GantryPosition,
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)

    colptr = D.colptr
    rowval = D.rowval
    nzval = D.nzval

    # Fill colptr
    fill_colptr!(D, pos, bixels, gantry, calc)

    # Preallocate arrays
    nprealloc = colptr[end]-1
    resize!(nzval, nprealloc)
    resize!(rowval, nprealloc)

    # Compute row and matrix values
    dose_kernel!(D, pos, bixels, gantry, surf, calc)

    D
end

function dose_fluence_matrix!(D::SparseMatrixCSC, pos, beamlets::AbstractVector{<:Beamlet},
    surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm)

    colptr = D.colptr
    rowval = D.rowval
    nzval = D.nzval

    # Fill colptr
    fill_colptr!(D, pos, beamlets, calc)

    # Preallocate arrays
    nprealloc = colptr[end]-1
    resize!(nzval, nprealloc)
    resize!(rowval, nprealloc)

    # Compute row and matrix values
    dose_kernel!(D, pos, beamlets, surf, calc)

    D
end
