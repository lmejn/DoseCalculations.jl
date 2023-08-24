#
#   DICOM Plan
#
# Functions for loading in DICOM Plan data. Currently only tested on VMAT SPARK
# data.
#

export load_dicom

#--- DICOM IO ----------------------------------------------------------------------------------------------------------

"""
    load_dicom(filename)

Load a DICOM RP file into a Vector{TreatmentField}.
"""
function load_dicom(filename)
    dcm = dcm_parse(filename)

    referenced_dose = load_ref_dose.(dcm[tag"FractionGroupSequence"][1].ReferencedBeamSequence)
    load_beam.(dcm[tag"BeamSequence"], referenced_dose)
end

"""
    load_beam(beam, total_meterset)

Load a beam from a control point sequence in a DICOM RP file.
"""
function load_beam(beam, total_meterset)

    SAD = beam[tag"SourceAxisDistance"]

    controlpoints = beam[tag"ControlPointSequence"]

    controlpoint = controlpoints[1]

    # 
    ncontrol = beam[tag"NumberOfControlPoints"]

    θb = deg2rad(controlpoint[tag"BeamLimitingDeviceAngle"])
    doserate = fill(controlpoint[tag"DoseRateSet"]/60., ncontrol) # Convert from MU/min to MU/s

    isocenter = fill(SVector(controlpoint[tag"IsocenterPosition"]...), ncontrol)

    # Jaws
    jaws_x = controlpoint[tag"BeamLimitingDevicePositionSequence"][1][tag"LeafJawPositions"]
    jaws_y = controlpoint[tag"BeamLimitingDevicePositionSequence"][2][tag"LeafJawPositions"]
    jaws = fill(Jaws(jaws_x, jaws_y), ncontrol)

    mlcy = beam[tag"BeamLimitingDeviceSequence"][3]["LeafPositionBoundaries"]
    nleaves = length(mlcy)-1

    T = typeof(θb)

    ϕg = zeros(ncontrol)
    meterset = zeros(ncontrol)
    sourcepositions = Vector{RotatingGantryPosition{T}}(undef, ncontrol)

    mlcs = Vector{MultiLeafCollimator{Matrix{T}, Vector{T}}}(undef, ncontrol)

    for i in eachindex(controlpoints)

        controlpoint = controlpoints[i]

        ϕg = fixangle(deg2rad(controlpoint[tag"GantryAngle"]))
        sourcepositions[i] = RotatingGantryPosition(ϕg, θb, SAD)

        meterset[i] = total_meterset*controlpoint[tag"CumulativeMetersetWeight"]

        mlcx = reshape(controlpoint[tag"BeamLimitingDevicePositionSequence"][end][tag"LeafJawPositions"], nleaves, 2)'
        mlcs[i] = MultiLeafCollimator(Array(mlcx), mlcy)
    end

    TreatmentField(mlcs, sourcepositions, meterset, doserate, isocenter)
end

"""
    load_ref_dose(beam)

Load a reference dose, used for calculating the meterset in a control point
sequence
"""
load_ref_dose(beam) = beam[tag"BeamMeterset"]
