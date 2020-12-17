struct SphereSignedDistanceFunction{T<:Real} <: AbstractSignedDistanceFunction
    radius::T
end

radius(sdf::SphereSignedDistanceFunction) = sdf.radius
(sdf::SphereSignedDistanceFunction{T1})(point::Point{3, T2}) where {T1<:Real, T2<:Real} = 
    norm(point) - T2(radius(sdf))

# Currently always pointing out 
normal(::SphereSignedDistanceFunction{T1}, point::Point{3, T2}) where {T1<:Real, T2<:Real} = 
    Vec3{T2}(normalize(point))