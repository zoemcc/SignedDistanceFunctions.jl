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

overshootfactor(bvsdf::BoundingVolumeSignedDistanceFunction) = bvsdf.overshootfactor
boundingSDF(bvsdf::BoundingVolumeSignedDistanceFunction) = bvsdf.boundingSDF
interiorSDF(bvsdf::BoundingVolumeSignedDistanceFunction) = bvsdf.interiorSDF


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

SDFs(union_sdf::UnionSignedDistanceFunction) = union_sdf.SDFs

function (union_sdf::UnionSignedDistanceFunction)(point::Point{3, T}) where {T<:Real}
    minimum(sdf(point) for sdf in SDFs(union_sdf))
end

function normal(union_sdf::UnionSignedDistanceFunction, point::Point{3, T}) where {T<:Real}
    min_sdf_index = argmin(sdf(point) for sdf in SDFs(union_sdf))
    normal(SDFs(union_sdf)[min_sdf_index], point)
end

struct IntersectSignedDistanceFunction{T<:Tuple{Vararg{AbstractSignedDistanceFunction}}} <: AbstractSignedDistanceFunction
    SDFs::T
end

SDFs(intersect_sdf::IntersectSignedDistanceFunction) = intersect_sdf.SDFs

function (intersect_sdf::IntersectSignedDistanceFunction)(point::Point{3, T}) where {T<:Real}
    maximum(sdf(point) for sdf in SDFs(intersect_sdf))
end

function normal(intersect_sdf::IntersectSignedDistanceFunction, point::Point{3, T}) where {T<:Real}
    max_sdf_index = argmax(sdf(point) for sdf in SDFs(intersect_sdf))
    normal(SDFs(intersect_sdf)[max_sdf_index], point)
end

struct NegatedSignedDistanceFunction{T<:AbstractSignedDistanceFunction} <: AbstractSignedDistanceFunction
    sdf::T
end

sdf(negated_sdf::NegatedSignedDistanceFunction) = negated_sdf.sdf

function (negated_sdf::NegatedSignedDistanceFunction)(point::Point{3, T}) where {T<:Real}
    -(sdf(negated_sdf)(point))
end

function normal(negated_sdf::NegatedSignedDistanceFunction, point::Point{3, T}) where {T<:Real}
    -(normal(sdf(negated_sdf), point))
end


struct TransformedSignedDistanceFunction{Transform<:Transformation, SDF<:AbstractSignedDistanceFunction} <: AbstractSignedDistanceFunction
    transform::Transform
    invtransform::Transform
    sdf::SDF
end

transform(transformed_sdf::TransformedSignedDistanceFunction) = transformed_sdf.transform
invtransform(transformed_sdf::TransformedSignedDistanceFunction) = transformed_sdf.invtransform
sdf(transformed_sdf::TransformedSignedDistanceFunction) = transformed_sdf.sdf
transform_sdf(transformation::Transformation, sdf::AbstractSignedDistanceFunction) = TransformedSignedDistanceFunction(transformation, inv(transformation), sdf)
function transform_sdf(transformation::Transformation, transformed_sdf::TransformedSignedDistanceFunction) 
    # collapse composed transformations
    composed_transformation = compose(transformation, transform(transformed_sdf))
    TransformedSignedDistanceFunction(composed_transformation, inv(composed_transformation), sdf(transformed_sdf))
end

(transformed_sdf::TransformedSignedDistanceFunction)(point::Point{3, T}) where {T<:Real} = sdf(transformed_sdf)(invtransform(transformed_sdf)(point))

function normal(transformed_sdf::TransformedSignedDistanceFunction, point::Point{3, T}) where {T<:Real}
    Vec{3, T}(normalize(LinearMap(transform_deriv(invtransform(transformed_sdf), invtransform(transformed_sdf)(point)))(invtransform(transformed_sdf)(point))))
end


