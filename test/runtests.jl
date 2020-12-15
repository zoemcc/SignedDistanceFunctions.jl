using SignedDistanceFunctions
using GeometryBasics
using CoordinateTransformations
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

    # translation transformation tests
    translation1 = Translation(1., 0., 0.)
    translated_sphere1 = SignedDistanceFunctions.transform_sdf(translation1, spheresdf1)
    translated_sphere2 = SignedDistanceFunctions.transform_sdf(translation1, spheresdf2)
    @test translated_sphere1(p1) ≈ 0.5
    @test translated_sphere1(p2) ≈ -0.2
    @test translated_sphere2(p1) ≈ 0.0
    @test translated_sphere2(p2) ≈ -0.7

    translation2 = Translation(-0.5, 0., 0.)
    translated2x_sphere1 = SignedDistanceFunctions.transform_sdf(translation2, translated_sphere1)
    translated2x_sphere2 = SignedDistanceFunctions.transform_sdf(translation2, translated_sphere2)
    @test translated2x_sphere1(p1) ≈ 0.0
    @test translated2x_sphere1(p2) ≈ 0.3
    @test translated2x_sphere2(p1) ≈ -0.5
    @test translated2x_sphere2(p2) ≈ -0.2





end
