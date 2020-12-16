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
include("sdfs.jl")

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

    TransformedSignedDistanceFunction,
    transform, invtransform, sdf, transform_sdf,

    SphereSignedDistanceFunction,
    radius
    


end
