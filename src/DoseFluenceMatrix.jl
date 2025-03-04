export dose_fluence_matrix, dose_fluence_matrix!

#--- Dose-Fluence Matrix-------------------------------------------------------

"""
    dose_fluence_matrix(T, pos, beamlets, surf, calc)

Compute a dose-fluence matrix from dose positions, beamlets, external surface and
dose calculation algorithm.

`T` is the matrix type. It currently supports:
- `Matrix`: Dense CPU
- `SparseMatrixCSC`: Sparse CPU
- `CuArray`: Dense GPU
e.g. `dose_fluence_matrix(SparseMatrixCSC, ...)` will create a sparse matrix
computed using the CPU.

See `dose_fluence_matrix!` for implementation.
"""
function dose_fluence_matrix(T::Type{<:AbstractMatrix}, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                             kwargs...)
    D = zeros(length(pos), length(beamlets))
    dose_fluence_matrix!(D, pos, beamlets, surf, calc; kwargs...)
end

function dose_fluence_matrix(::Type{CuArray}, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                             kwargs...)
    D = CUDA.zeros(length(pos), length(beamlets))
    dose_fluence_matrix!(D, pos, beamlets, surf, calc; kwargs...)
end

function dose_fluence_matrix(::Type{SparseMatrixCSC}, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                             surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                             kwargs...)
    D = spzeros(length(pos), length(beamlets))
    dose_fluence_matrix!(D, pos, beamlets, surf, calc; kwargs...)
end

"""
    dose_fluence_matrix!(D<:AbstractMatrix, pos, beamlets, surf, calc)

Compute a dose-fluence matrix from dose positions, beamlets, external surface
and dose calculation algorithm.

Requires the `point_dose` method to be defined for the given dose calculation
algorithm (`calc`). `point_dose` computes the dose at a given position from a
given beamlet.

It stores result in `D<:AbstractMatrix`. Currently, the following matrix types
are supported:
- `D::Matrix`
- `D::SparseMatrixCSC`
- `D::CuArray`

"""
dose_fluence_matrix!

function dose_fluence_matrix!(D::SparseMatrixCSC, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                              maxradius=100.)
    _assert_size(D, pos, beamlets)

    colptr = D.colptr
    rowval = D.rowval
    nzval = D.nzval

    # Fill colptr
    fill_colptr!(D, pos, beamlets, maxradius)

    # Preallocate arrays
    nprealloc = colptr[end]-1
    resize!(nzval, nprealloc)
    resize!(rowval, nprealloc)

    # Fill rowval
    fill_rowval!(D, pos, beamlets, maxradius)

    # Compute row and matrix values
    dose_kernel!(D, pos, beamlets, surf, calc)

    D
end

function dose_fluence_matrix!(D::AbstractArray, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                              maxradius=100.)
    _assert_size(D, pos, beamlets)
    @batch per=thread for j in eachindex(beamlets), i in eachindex(pos)
        @inbounds D[i, j] = point_dose(pos[i], beamlets[j], surf, calc, maxradius)
    end
    D
end

function dose_fluence_matrix!(D::CuArray, pos, beamlets::AbstractArray{<:AbstractBeamlet},
                              surf::AbstractExternalSurface, calc::AbstractDoseAlgorithm;
                              maxradius=100.)
    _assert_size(D, pos, beamlets)
    D .= point_dose.(vec(pos), permutedims(vec(beamlets)), Ref(surf), Ref(calc), maxradius)
end

"""
    _assert_size(D, pos, beamlets)

Ensures size of D is correct
"""
_assert_size(D, pos, beamlets) = @assert size(D) == (length(pos), length(beamlets))

function point_dose(pos::SVector{3, T}, beamlet, surf, calc, maxradius) where T<:Number
    src = source_position(beamlet)
    a = direction(beamlet)
    SAD = source_axis_distance(beamlet)  

    !kernel_size(pos-src, a, maxradius/SAD) && return zero(T)
    point_dose(pos, beamlet, surf, calc)
end

#--- Dose-Fluence Matrix Computations -----------------------------------------

kernel_size(r::SVector{3}, a::SVector{3}, maxradius) = sum(r.^2) < sum(r.*a)^2*(1+maxradius^2)

function fill_colptr!(D::SparseMatrixCSC, pos, beamlets, maxradius)
    colptr = D.colptr

    colptr[1] = 1
    @batch per=thread for j in eachindex(beamlets)
        beamlet = beamlets[j]

        src = source_position(beamlet)
        a = direction(beamlet)
        SAD = source_axis_distance(beamlet)    

        n = 0
        for i in eachindex(pos)
            n += kernel_size(pos[i]-src, a, maxradius/SAD)
        end
        colptr[j+1] = n
    end
    cumsum!(colptr, colptr)
end

function fill_rowval!(D::SparseMatrixCSC, pos, beamlets, maxradius)

    colptr = D.colptr
    rowval = D.rowval

    @batch per=thread for j in eachindex(beamlets) #
        beamlet = beamlets[j]

        # Create views of rowval and nzval
        ptr = colptr[j]:(colptr[j+1]-1)
        I = @view rowval[ptr]

        s = source_position(beamlet)
        a = direction(beamlet)
        SAD = source_axis_distance(beamlet)  

        n = 0
        for i in eachindex(pos)
            p = pos[i]
            if kernel_size(p-s, a, maxradius/SAD)
                n += 1
                I[n] = i
            end
        end
    end

end

"""
    dose_kernel!(rowval, nzval, pos::AbstractVector{T}, bixels, surf, calc)

Compute the fluence kernel for a given position.

Designed to be used with a dose-fluence matrix of type SparseMatrixCSC. Stores
the row value in `rowval`, and dose value in `nzval`.
"""
function dose_kernel!(D::SparseMatrixCSC, pos, beamlets, surf, calc)

    colptr = D.colptr
    rowval = D.rowval
    nzval = D.nzval

    jprev = 1
    @batch per=thread for n in eachindex(rowval, nzval)
        i = rowval[n]
        j = _sequential_searchsortedlast(colptr, n, jprev)
        nzval[n] = point_dose(pos[i], beamlets[j], surf, calc)
        jprev = j
    end

end

"""
    _sequential_searchsortedlast(a, x, jstart)

A `searchsortedlast` where the index is near and greater than previous index

This assumes that the next index is larger than the starting index `jstart`. It first
checks that `jstart` already satisfies `a[j]<=x<a[j+1]`. If not, it searches for
`j>jstart`, returning `length(a)` if `a[end]<x`.
"""
function _sequential_searchsortedlast(a, x, jstart)
    @inbounds a[jstart]<=x<a[jstart+1] && return jstart
    for j = jstart+1:length(a)-1
        @inbounds a[j]<=x<a[j+1] && return j
    end
    length(a)
end
