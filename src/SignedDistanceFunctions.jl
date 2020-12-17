module SignedDistanceFunctions

using GeometryBasics
using Rotations
using CoordinateTransformations
using StaticArrays
using LinearAlgebra

using FiniteDifferences
using ForwardDiff
using Zygote

include("types.jl")
include("shapes.jl")

export 
    AbstractSignedDistanceFunction,

    normal,
    normalforwarddiff,
    normalzygote,
    normalfinitediff,

    BoundingVolumeSignedDistanceFunction, 
    overshootfactor, boundingSDF, interiorSDF,

    UnionSignedDistanceFunction,
    SDFs,

    IntersectSignedDistanceFunction,

    NegatedSignedDistanceFunction,
    sdf,

    TransformedSignedDistanceFunction,
    transform, invtransform, transform_sdf,

    SphereSignedDistanceFunction,
    radius
    


end
