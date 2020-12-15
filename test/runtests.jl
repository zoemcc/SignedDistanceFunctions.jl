using SignedDistanceFunctions
using GeometryBasics
using Test

@testset "SignedDistanceFunctions.jl" begin
    T = Float64
    spheresdf1 = SignedDistanceFunctions.SphereSignedDistanceFunction(0.5)
    spheresdf2 = SignedDistanceFunctions.SphereSignedDistanceFunction(1.0)
    unionsdf = SignedDistanceFunctions.UnionSignedDistanceFunction((spheresdf1, spheresdf2))

    p1 = zero(Point{3, T})
    @test spheresdf1(p1) ≈ -0.5
    @test spheresdf2(p1) ≈ -1.0
    @test unionsdf(p1) ≈ -1.0

    bvsdf = SignedDistanceFunctions.BoundingVolumeSignedDistanceFunction(0.25, spheresdf2, spheresdf1)

    p2 = Point{3, T}(1.3, 0., 0.)
    p3 = Point{3, T}(1.2, 0., 0.)
    @test bvsdf(p1) ≈ -0.5
    @test bvsdf(p2) ≈ 0.3
    @test bvsdf(p3) ≈ 0.7




end
