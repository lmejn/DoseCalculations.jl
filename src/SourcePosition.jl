export RotatingGantryPosition

#--- Source Position ----------------------------------------------------------

abstract type AbstractSourcePosition{T} <: AbstractVector{T} end

Base.size(source::AbstractSourcePosition) = (3,)

#--- Rotating Gantry Position -------------------------------------------------

"""
    RotatingGantryPosition{T}

The position/rotation of the source and beam-limiting device.

Stores:
- `gantry_angle`: As defined by IEC
- `collimator_angle`: As defined by IEC (beam limiting device angle)
- `source_axis_distance`: Distance between source and isocenter
- `central_beam_axis`: Unit vector pointing from the isocenter to the source in IEC Fixed coordinates.
"""
struct RotatingGantryPosition{T} <: AbstractSourcePosition{T}
    gantry_angle::T
    collimator_angle::T
    source_axis_distance::T
    central_beam_axis::SVector{3, T}

    function RotatingGantryPosition(ϕg, θb, SAD)
        ϕg, θb, SAD = promote(ϕg, θb, SAD)
        T = typeof(ϕg)
        ax = SVector(sin(ϕg), zero(T), cos(ϕg))
        new{T}(ϕg, θb, SAD, ax)
    end
end

function Base.show(io::IO, source::RotatingGantryPosition)
    @printf io "RotatingGantryPosition(ϕg=%0.1f°, θb=%0.1f°, SAD=%0.1f)" rad2deg(getϕg(source)) rad2deg(getθb(source)) getSAD(source)
end

Base.getindex(source::RotatingGantryPosition, i::Int) = getSAD(source)*beamaxis(source, i)

getϕg(source::RotatingGantryPosition) = source.gantry_angle
getθb(source::RotatingGantryPosition) = source.collimator_angle
getSAD(source::RotatingGantryPosition) = source.source_axis_distance
beamaxis(source::RotatingGantryPosition) = source.central_beam_axis
beamaxis(source::RotatingGantryPosition, i) = source.central_beam_axis[i]

getposition(source::RotatingGantryPosition) = getSAD(source)*beamaxis(source)

fixed_to_bld(source::RotatingGantryPosition) = fixed_to_bld(getϕg(source), getθb(source), getSAD(source))
bld_to_fixed(source::RotatingGantryPosition) = bld_to_fixed(getϕg(source), getθb(source), getSAD(source))
