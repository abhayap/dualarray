using Test

include("../src/modalcoeff.jl")

@testset "modal coefficient" begin
    @testset "open sphere" begin
        @test bnopen(0, 0.1) ≈ 12.545437132817783 + 0.0im atol=1e-14
        @test bnopen(1, 0.1) ≈ 0.0 + 0.41846029103011073im atol=1e-14
        @test bnopen(2, 0.1) ≈ -0.00837159808553246 + 0.0im atol=1e-14
    end
    @testset "rigid sphere" begin
        @test bnrigid(0, 0.1, 0.1) ≈ 12.504005420834602 + 0.004143171198317932im atol=1e-14
        @test bnrigid(1, 0.1, 0.1) ≈ 0.00010440317272607134 + 0.6283106682095139im atol=1e-14
        @test bnrigid(2, 0.1, 0.1) ≈ -0.013954900578375413 + 2.0632970657325567e-9im atol=1e-14
        @test bnrigid(0, 0.0, 0) ≈ 4π atol=1e-14
        @test bnrigid(1, 0.0, 0) ≈ 0 atol=1e-14
        @test bnrigid.(0:3, 1.3, 1.3) ≈ [7.101281585795888 + 2.876763204917208im,
                                         0.945455469319844 + 6.166942016106375im,
                                         -2.173494107438295 + 0.085477251197533im,
                                         -0.000999409593938 - 0.425087990099182im] atol=1e-14
    end
    @testset "dual sphere" begin
        @test bnrigid(0, 3, 0.4) ≈ 0.6717145683184501 + 0.013095809333291006im atol=1e-14
        @test bnrigid(1, 3, 0.4) ≈ 0.04390812204821476 + 4.335454183794927im atol=1e-14
        @test bnrigid(2, 3, 0.4) ≈ -3.7532826771091092 + 0.0005517351494948208im atol=1e-14
    end
    @testset "cardioid sphere" begin
        @test bncardioid(0, 0.1) ≈ 12.545437132817783 + 0.41846029103011073im atol=1e-14
        @test bncardioid(1, 0.1) ≈ 4.176231312215574 + 0.41846029103011073im atol=1e-14
        @test bncardioid(2, 0.1) ≈ -0.00837159808553246 + 0.16731234846413676im atol=1e-14
    end
end
