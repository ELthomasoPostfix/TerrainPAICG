# Utils that do NOT require a noise map argument.
@testset "Utils" begin
    @testset "Cartesian Product Ordering" begin
        lists = [[2,3],[1,2],[3,4]]
        cartprod = Iterators.product(lists...)
        cartprod = reduce(vcat, cartprod) # Reduce the matrix to a vector.
        # If the product input is sorted, then the product output may be sorted.
        @test all(issorted(lst) for lst in lists)
        # Test that `Iterators.product` outputs tuples in an expected ordering.
        @test cartprod == [
            (2, 1, 3), (3, 1, 3), (2, 2, 3), (3, 2, 3),
            (2, 1, 4), (3, 1, 4), (2, 2, 4), (3, 2, 4)
        ]
        # That expected ordering must be a sorted ordering.
        # Given tuples (a, b, c) we expect them to be sorted by priority c > b > a.
        # So, if two tuples have the same c, then sort them by b. If they have the
        # same c and b, then sort them by a.
        @test issorted(cartprod, by=reverse)

        # Ensure the utils function produces the same result as tested above.
        @test cartprod == GeneralNoise.corners((2, 1, 3))
    end

    @testset "Linear Interpolation" begin
        a = -2.0
        b = +5.0
        # Test out points.
        @test_throws AssertionError GeneralNoise.lerp(0.0 - 1e-10, a, b)
        @test_throws AssertionError GeneralNoise.lerp(1.0 + 1e-10, a, b)
        # Test boundary points.
        @test GeneralNoise.lerp(0.0, a, b) == a
        @test GeneralNoise.lerp(1.0, a, b) == b
        # Test in points.
        @test isapprox(GeneralNoise.lerp(0.1, a, b), -1.3, atol=1e-4)
        @test isapprox(GeneralNoise.lerp(0.3, a, b), +0.1, atol=1e-4)
        @test isapprox(GeneralNoise.lerp(0.5, a, b), +1.5, atol=1e-4)
        @test isapprox(GeneralNoise.lerp(0.7, a, b), +2.9, atol=1e-4)
        @test isapprox(GeneralNoise.lerp(0.9, a, b), +4.3, atol=1e-4)
        # Test the range output range thoroughly, it must be [a, b].
        @test all(f -> a <= f <= b , [ GeneralNoise.lerp(x, a, b) for x in 0:0.01:1 ])
    end

    @testset "Smoothstep Smoothing" begin
        # Test out points.
        @test_throws AssertionError GeneralNoise.smoothstep(0.0 - 1e-10)
        @test_throws AssertionError GeneralNoise.smoothstep(1.0 + 1e-10)
        # Test boundary points.
        @test GeneralNoise.smoothstep(0.0) == 0.0
        @test GeneralNoise.smoothstep(1.0) == 1.0
        # Test in points.
        @test isapprox(GeneralNoise.smoothstep(0.1), 0.028, atol=1e-4)
        @test isapprox(GeneralNoise.smoothstep(0.3), 0.216, atol=1e-4)
        # Require exact equality for x=0.5. We may rely on this in other tests.
        @test GeneralNoise.smoothstep(0.5) == 0.5
        @test isapprox(GeneralNoise.smoothstep(0.7), 0.784, atol=1e-4)
        @test isapprox(GeneralNoise.smoothstep(0.9), 0.972, atol=1e-4)
        # Test the range output range thoroughly, it must be [0, 1].
        @test all(f -> 0.0 <= f <= 1.0 , [ GeneralNoise.smoothstep(x) for x in 0:0.01:1 ])
    end

    @testset "x = 2^n" begin
        @test ! ispow2(0)
        @test   ispow2(1) # 2^0
        @test   ispow2(2) # 2^1
        @test ! ispow2(3)
        @test   ispow2(4) # 2^2
        @test ! ispow2(5)
        @test ! ispow2(6)
        @test ! ispow2(7)
        @test   ispow2(8) # 2^3
        @test ! ispow2(9)
        @test ! ispow2(12)
        @test ! ispow2(15)
        @test   ispow2(16) # 2^4
        @test ! ispow2(17)
        @test ! ispow2(463)
        @test ! ispow2(1023)
        @test   ispow2(1024) # 2^10
    end
end

# Test constructor functions.
@testset "Constructors" begin
    # Test against invalid map dimensions..
    @test_throws AssertionError PerlinNoiseMap{0,Float64}(())
    @test_throws AssertionError PerlinNoiseMap{1,Float64}((1,))
    @test_throws AssertionError PerlinNoiseMap{Float64}(())
    @test_throws AssertionError PerlinNoiseMap{Float64}((1,))

    # The happy day constructor scenario.
    @testset "Happy Day" begin
        for dim in 2:4
            # Try creating a map with as few gridcells as possible.
            celldims = Tuple(fill(1, dim))  # (1, ..., 1)
            graddims = celldims .+ 1        # (2, ..., 2)
            mapInner = PerlinNoiseMap{dim,Float64}(celldims)
            @test size(mapInner.gradients) == graddims
            mapOuter = PerlinNoiseMap{Float64}(celldims)
            @test size(mapOuter.gradients) == graddims

            # Try creating a map with increasingly many gridcells.
            celldims = Tuple(fill(dim, dim))  # (d,   ..., d  )
            graddims = celldims .+ 1          # (d+1, ..., d+1)
            mapInner = PerlinNoiseMap{dim,Float64}(celldims)
            @test size(mapInner.gradients) == graddims
            mapOuter = PerlinNoiseMap{Float64}(celldims)
            @test size(mapOuter.gradients) == graddims
        end
    end

    # A dimension of the map can not be <= 0 gridcells wide.
    @test_throws AssertionError PerlinNoiseMap{Float64}((0,  1))
    @test_throws AssertionError PerlinNoiseMap{Float64}((1,  0))
    @test_throws AssertionError PerlinNoiseMap{Float64}((1, -1))
    @test_throws AssertionError PerlinNoiseMap{Float64}((2, -100, 3))
end

# Utils that DO require a noise map argument.
@testset "Map Utils" begin
    map = PerlinNoiseMap{Float64}((1, 4, 3))

    @test GeneralNoise.dims(map) == (2, 5, 4)
    @test GeneralNoise.gridmax(map) == (1, 4, 3)

    # Every map contains at least one gridcell, so (0, 0, 0) is always in bounds.
    @test GeneralNoise.inbounds(map, (0, 0, 0))
    # Test the largest possible in bounds cell.
    @test GeneralNoise.inbounds(map, (0, 3, 2))
    # Test different ways to go out of bounds of the noise map.
    @test !GeneralNoise.inbounds(map, (-1, 0, 0))
    @test !GeneralNoise.inbounds(map, (1, 0, 0))
    @test !GeneralNoise.inbounds(map, (0, 4, 0))
    @test !GeneralNoise.inbounds(map, (0, 0, 3))

    # Try to get the gridcell indices matching sample points.
    @test GeneralNoise.gridcell(map, (0.0, 0.0, 0.0)) == (0, 0, 0)
    @test GeneralNoise.gridcell(map, (0.9, 3.9, 2.9)) == (0, 3, 2)
    @test GeneralNoise.gridcell(map, (0.9999, 3.9999, 2.9999)) == (0, 3, 2)
    @test GeneralNoise.gridcell(map, (1.0, 4.0, 3.0)) == (0, 3, 2)
    @test_throws AssertionError GeneralNoise.gridcell(map, (1.1, 4.1, 3.1))
    @test_throws AssertionError GeneralNoise.gridcell(map, (1.0001, 4.0001, 3.0001))
    @test_throws AssertionError GeneralNoise.gridcell(map, (0.0, -0.1, 0.0))
end

# Verify that the bounds of the noise range are as expected.
# I.e. we know under which circumstances the highest and lowest possible
# perlin noise values manifest, and for which sample point they do.
# Given a n-dimensional noise map, all noise perlin values are naturally
# bounded in the range [-sqrt(n)/2, sqrt(n)/2].
@testset "Noise Range Bounds" begin
    @testset "1x1 map max noise" begin
        # Create a 2D grid with one single, square gridcell.
        map = PerlinNoiseMap{Float64}((1, 1))

        # All gradient vectors point directly TOWARDS the exact gridcell center.
        map.gradients[1,1] = Tuple(normalize([+1, +1])) # bottom left  corner
        map.gradients[2,1] = Tuple(normalize([-1, +1])) # bottom right corner
        map.gradients[1,2] = Tuple(normalize([+1, -1])) # top    left  corner
        map.gradients[2,2] = Tuple(normalize([-1, -1])) # top    right corner
        result = noise(map, (0.5, 0.5))
        actual = + (sqrt(2) / 2)
        # This gradient + sample point configuration produces the HIGHEST possible noise.
        @test isapprox(result, actual, atol=1e-10)
    end

    @testset "1x1 map min noise" begin
        # Create a 2D grid with one single, square gridcell.
        map = PerlinNoiseMap{Float64}((1, 1))

        # All gradient vectors point directly AWAY FROM the exact gridcell center.
        map.gradients[1,1] = Tuple(normalize([-1, -1])) # bottom left  corner
        map.gradients[2,1] = Tuple(normalize([+1, -1])) # bottom right corner
        map.gradients[1,2] = Tuple(normalize([-1, +1])) # top    left  corner
        map.gradients[2,2] = Tuple(normalize([+1, +1])) # top    right corner
        result = noise(map, (0.5, 0.5))
        actual = - (sqrt(2) / 2)
        # This gradient + sample point configuration produces the LOWEST possible noise.
        @test isapprox(result, actual, atol=1e-10)
    end

    @testset "1x1x1 map max noise" begin
        # Create a 3D grid with one single, cubic gridcell.
        map = PerlinNoiseMap{Float64}((1, 1, 1))

        # All gradient vectors point directly TOWARDS the exact gridcell center.
        map.gradients[1,1,1] = Tuple(normalize([+1, +1, +1]))
        map.gradients[2,1,1] = Tuple(normalize([-1, +1, +1]))
        map.gradients[1,2,1] = Tuple(normalize([+1, -1, +1]))
        map.gradients[2,2,1] = Tuple(normalize([-1, -1, +1]))
        map.gradients[1,1,2] = Tuple(normalize([+1, +1, -1]))
        map.gradients[2,1,2] = Tuple(normalize([-1, +1, -1]))
        map.gradients[1,2,2] = Tuple(normalize([+1, -1, -1]))
        map.gradients[2,2,2] = Tuple(normalize([-1, -1, -1]))
        result = noise(map, (0.5, 0.5, 0.5))
        actual = + (sqrt(3) / 2)
        # This gradient + sample point configuration produces the HIGHEST possible noise.
        @test isapprox(result, actual, atol=1e-10)
    end

    @testset "1x1x1 map min noise" begin
        # Create a 3D grid with one single, cubic gridcell.
        map = PerlinNoiseMap{Float64}((1, 1, 1))

        # All gradient vectors point directly AWAY FROM the exact gridcell center.
        map.gradients[1,1,1] = Tuple(normalize([-1, -1, -1]))
        map.gradients[2,1,1] = Tuple(normalize([+1, -1, -1]))
        map.gradients[1,2,1] = Tuple(normalize([-1, +1, -1]))
        map.gradients[2,2,1] = Tuple(normalize([+1, +1, -1]))
        map.gradients[1,1,2] = Tuple(normalize([-1, -1, +1]))
        map.gradients[2,1,2] = Tuple(normalize([+1, -1, +1]))
        map.gradients[1,2,2] = Tuple(normalize([-1, +1, +1]))
        map.gradients[2,2,2] = Tuple(normalize([+1, +1, +1]))
        result = noise(map, (0.5, 0.5, 0.5))
        actual = - (sqrt(3) / 2)
        # This gradient + sample point configuration produces the LOWEST possible noise.
        @test isapprox(result, actual, atol=1e-10)
    end
end