import Base.getindex, Base.lastindex, Base.length, Base.iterate, Base.show

export getmlc, getjaws, getϕg, getθb, getmeterset, getdoserate, getisocenter, getSAD, getmetersettype
export resample, getgantry, getΔMU
export VMATField

fixangle(angle) = angle > π ? angle - 2π : angle

#--- Abstract Types ----------------------------------------------------------------------------------------------------

abstract type AbstractControlPoint end
const AbstractTreatmentField = AbstractVector{AbstractControlPoint}
const AbstractTreatmentPlan = AbstractVector{AbstractTreatmentField}

#--- ControlPoint ------------------------------------------------------------------------------------------------------

"""
    ControlPoint

Elements of a TreatmentField
"""
struct ControlPoint{T<:AbstractFloat, TBLD<:AbstractBeamLimitingDevice,
    TSourcePosition<:AbstractSourcePosition} <: AbstractControlPoint
    # Beam Limiting Devices
    beam_limiting_device::TBLD

    # Source Position
    source_position::TSourcePosition

    # Machine Parameters
    ΔMU::T # Meterset (MU)
    doserate::T # Dose Rate (MU/s)

    isocenter::SVector{3, T}    # Isocenter Position
end

function Base.show(io::IO, point::ControlPoint)
    @printf io "ΔMU=%0.1f" getΔMU(point)
    println(io, ", ", point.source_position)
    println(io, point.beam_limiting_device)
end

# Accessors
getBLD(point::ControlPoint) = point.beam_limiting_device
getsourceposition(point::ControlPoint) = point.source_position
getΔMU(point::ControlPoint) = point.ΔMU
getdoserate(point::ControlPoint) = point.doserate
getisocenter(point::ControlPoint) = point.isocenter

fixed_to_bld(point::ControlPoint) = point |> getsourceposition |> fixed_to_bld
getΔt(point::ControlPoint) = getΔMU(point)/getdoserate(point)

#--- TreatmentField ----------------------------------------------------------------------------------------------------

struct TreatmentField{TBLD<:AbstractBeamLimitingDevice, TSrcPos<:AbstractSourcePosition,
    T<:Real} <: AbstractTreatmentField
    n::Int

    beamlimitingdevices::Vector{TBLD}
    sourcepositions::Vector{TSrcPos}
    
    # Machine Parameters
    meterset::Vector{T} # Cumulative Meterset (MU)
    doserate::Vector{T} # Dose Rate (MU/s)

    isocenter::Vector{SVector{3, T}}    # Isocenter Position

    function TreatmentField(beamlimitingdevices, sourcepositions, meterset, doserate, isocenter)
        n = length(beamlimitingdevices)
        @assert length(sourcepositions) == n "All input vectors must be of the same length"
        @assert length(doserate) == n "All input vectors must be of the same length"
        @assert length(meterset) == n "All input vectors must be of the same length"
        @assert length(isocenter) == n "All input vectors must be of the same length"

        TBLD = eltype(beamlimitingdevices)
        TSrcPos = eltype(sourcepositions)
        T = eltype(meterset)
        new{TBLD, TSrcPos, T}(n, beamlimitingdevices, sourcepositions, meterset, doserate, isocenter)
    end
end

Base.size(field::TreatmentField) = (field.n,)

function Base.getindex(field::TreatmentField, i::Int)
    n = length(field)
    MU = field.meterset

    ΔMU = 0.5*(MU[min(n, i+1)]-MU[max(1, i-1)])

    ControlPoint(
        field.beamlimitingdevices[i],
        field.sourcepositions[i],
        ΔMU,
        field.doserate[i],
        field.isocenter[i]
        )
end

#--- Resampling --------------------------------------------------------------------------------------------------------

"""
    resample(field, ηₛ::AbstractVector{T}; by=:time)

Resample a treatment field onto new times or meterset values.

Can resample either by time (`by=:time`) or MU (`by=:MU`, default).
"""
function resample(field::AbstractTreatmentField, ηₛ::AbstractVector{T}; by=:MU) where T<:AbstractFloat

    if(by==:time)
        η = field.meterset./field.dose_rate
    elseif(by==:MU)
        η = field.meterset
    end
    Δη = @. η[2:end]-η[1:end-1]

    ncontrol = length(ηₛ)
    ϕg = zeros(ncontrol)
    meterset = zeros(ncontrol)

    mlc = MultiLeafCollimatorSequence(ncontrol, getedges(field.mlc))

    for i in eachindex(ηₛ)
        j = min(length(Δη), searchsortedlast(η, ηₛ[i]))
        α = (ηₛ[i] - η[j])/Δη[j]

        # Interpolate Gantry Angle
        ϕg[i] = interp(field.gantry_angle[j], field.gantry_angle[j+1], α)
        meterset[i] = interp(field.meterset[j], field.meterset[j+1], α)

        # Interpolate MLC x
        mlc.positions[:, :, i] .= interp(field.mlc.positions[:,:,j], field.mlc.positions[:,:,j+1], α)        

    end
    VMATField(mlc, field.jaws,
              ϕg, field.collimator_angle, field.source_axis_distance,
              meterset, field.dose_rate,
              field.isocenter)
end

"""
    resample(field, Δη::T; by=:MU)

Resample at uniform steps `Δη`, from start to finish.

See `resample(field, ηₛ::AbstractVector{T}; by=:MU)` for details.
"""
function resample(field, Δη::T; by=:MU, include_end=true) where T<:Number

    if(by==:time)
        η_start = field.meterset[1]/field.dose_rate
        η_end = field.meterset[end]/field.dose_rate
    elseif(by==:MU)
        η_start = field.meterset[1]
        η_end = field.meterset[end]
    end

    η = collect(η_start:Δη:η_end)
    if(include_end && η[end] != η_end)
        push!(η, η_end)
    end
    resample(field, η; by=by)
end
