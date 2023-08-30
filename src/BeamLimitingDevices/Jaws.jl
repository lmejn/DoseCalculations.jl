
#--- Jaws --------------------------------------------------------------------------------------------------------------

abstract type AbstractJaws <: AbstractBeamLimitingDevice end

""" 
    Jaws

`Jaws` stores the x and y positions of the jaws.

The x/y positions of the jaws can be accessed through the `getx`/`gety` methods.

The usual constructor directly takes the jaw x and y position vectors, but a
single fieldsize can also be specified.
"""
struct Jaws{T} <: AbstractJaws
    x::SVector{2, T}
    y::SVector{2, T}

    # Constructors
    function Jaws(x::AbstractVector{T}, y::AbstractVector{T}) where T<:AbstractFloat
        new{T}(SVector{2, T}(x), SVector{2, T}(y))
    end
    function Jaws(x1::T, x2::T, y1::T, y2::T) where T<:AbstractFloat
        new{T}(SVector(x1, x2), SVector(y1, y2))
    end
end

"""
    Jaws(fieldsize::T)

Create jaws with a square field of length `fieldsize`.
"""
Jaws(fieldsize::T) where T<:AbstractFloat = Jaws(-0.5*fieldsize, 0.5*fieldsize, -0.5*fieldsize, 0.5*fieldsize)

getx(jaws::Jaws) = jaws.x
gety(jaws::Jaws) = jaws.y
getpositions(jaws::Jaws) = getx(jaws), gety(jaws)

getarea(jaws::Jaws) = (jaws.x[2]-jaws.x[1])*(jaws.y[2]-jaws.y[1])

import Base.intersect
function Base.intersect(rect1::Jaws, rect2::Jaws)
    x1, y1 = getpositions(rect1)
    x2, y2 = getpositions(rect2)

    Jaws(max(x1[1], x2[1]), min(x1[2], x2[2]), max(y1[1], y2[1]), min(y1[2], y2[2]))
end
