import Base.+, Base.-, Base.getindex, Base.lastindex, Base.size, Base.length
import Base.(==), Base.view

export Jaws, MultiLeafCollimator, getx, gety, inaperture, centerposition, edgeposition, extract_subset

"""
    Abstract type for Beam Limiting Devices are based on

`AbstractBeamLimitingDevice` contains the types for MultiLeafCollimators, Jaws, etc.
"""
abstract type AbstractBeamLimitingDevice end

#--- Implementations ---------------------------------------------------------------------------------------------------

include("Jaws.jl")
include("MultiLeafCollimator.jl")
include("MultiLeafCollimatorSequence.jl")
