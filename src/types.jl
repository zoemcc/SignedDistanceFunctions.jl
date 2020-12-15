abstract type AbstractSignedDistanceFunction end

(sdf::AbstractSignedDistanceFunction)(::Point{3, T}) where {T<:Real} = error("Signed Distance Function for $(typeof(sdf)) has not been defined.")
normal(sdf::AbstractSignedDistanceFunction, ::Point{3, T}) where {T<:Real} = error("Normal of Signed Distance Function for $(typeof(sdf)) has not been defined.")
normalzygote(sdf::AbstractSignedDistanceFunction, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(Zygote.gradient(sdf, point)[1]))
normalfinitediff(sdf::AbstractSignedDistanceFunction, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(FiniteDifferences.grad(central_fdm(5, 1), sdf, point)[1]))
normalforwarddiff(sdf::AbstractSignedDistanceFunction, point::Point{3, T}) where {T<:Real} = 
    Vec3{T}(normalize(ForwardDiff.gradient(sdf, point)))



struct BoundingVolumeSignedDistanceFunction{T<:Real, BoundingSDF<:AbstractSignedDistanceFunction,
                                         InteriorSDF<:AbstractSignedDistanceFunction} <:
                                         AbstractSignedDistanceFunction
    overshootfactor::T
    boundingSDF::BoundingSDF
    interiorSDF::InteriorSDF
end

overshootfactor(sdf::BoundingVolumeSignedDistanceFunction) = sdf.overshootfactor
boundingSDF(sdf::BoundingVolumeSignedDistanceFunction) = sdf.boundingSDF
interiorSDF(sdf::BoundingVolumeSignedDistanceFunction) = sdf.interiorSDF


function (bvsdf::BoundingVolumeSignedDistanceFunction)(point::Point{3, T}) where {T<:Real} 
    bounding_signed_distance = boundingSDF(bvsdf)(point)
    if bounding_signed_distance > T(overshootfactor(bvsdf))
        bounding_signed_distance
    else
        interiorSDF(bvsdf)(point)
    end
end

function normal(bvsdf::BoundingVolumeSignedDistanceFunction, point::Point{3, T}) where {T<:Real} 
    bounding_signed_distance = boundingSDF(bvsdf)(point)
    if bounding_signed_distance > T(overshootfactor(bvsdf))
        normal(boundingSDF(bvsdf), point)
    else
        normal(interiorSDF(bvsdf), point)
    end
end

struct UnionSignedDistanceFunction{T<:Tuple{Vararg{AbstractSignedDistanceFunction}}} <: AbstractSignedDistanceFunction
    SDFs::T
end

SDFs(sdf::UnionSignedDistanceFunction) = sdf.SDFs

function (union_sdf::UnionSignedDistanceFunction)(point::Point{3, T}) where {T<:Real}
    minimum(sdf(point) for sdf in SDFs(union_sdf))
end

function normal(union_sdf::UnionSignedDistanceFunction, point::Point{3, T}) where {T<:Real}
    min_sdf_index = argmin(sdf(point) for sdf in SDFs(union_sdf))
    normal(SDFs(union_sdf)[min_sdf_index], point)
end

