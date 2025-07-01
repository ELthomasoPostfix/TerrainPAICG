"""An N-dimensional lattice to facilitate sampling N-dimensional Perlin noise.

The map represents a fixed-size N-dimensional lattice of gradient vectors.
That is, at creation you specify how many gridcells each dimension contains,
and then gradient vectors are sampled immediately for all possible gridpoints.
For example, a 3x4x2 grid has three dimensions and in total 3*4*2=24 gridpoints,
so 24 gradient vectors are sampled.
"""
struct PerlinNoiseMap{N, T<:Number}
    """The N-dimensional gradient vector map."""
    gradients :: Array{NTuple{N, T}, N}

    """PerlinNoiseMap inner constructor.

    @param[in] cells The number of grid cells (hypercubes) in each dimension.
    """
    function PerlinNoiseMap{N, T}(cells :: NTuple{N, I}) where {N, T<:Number, I<:Integer}

        @assert N > 0 "The noise map must be at least 1-dimensional."
        @assert all( cells .> 0 ) "Require 1 or more gridcells in each dimension."

        # Each corner of a grid cell needs an associated gradient vector.
        # There is one more corner/gradiant in each dimension than gridcells.
        gradients = Array{NTuple{N, T}, N}(undef, cells .+ 1)
        # Generate unit gradient vectors for each gridpoint.
        for idx in eachindex(gradients)
            # A 1D gradient varies in [-1, 1], it is not normalized to 1 or -1.
            if N == 1
                gradients[idx] = Tuple(rand(-1:0.01:1, N))
            # Normalization is strictly required for vectors of > 1 dimensions.
            else
                gradients[idx] = Tuple(normalize(rand(-1:0.01:1, N)))
            end
        end

        return new{N, T}(gradients)
    end

end

"""PerlinNoiseMap outer constructor."""
function PerlinNoiseMap{T}(cells :: NTuple{N, I}) where {N, T<:Number, I}
    return PerlinNoiseMap{N, T}(cells)
end

"""Determine the number of gradients in each dimension."""
dims(map::PerlinNoiseMap{N, T}) where {N, T} = size(map.gradients)

"""Determine the dimensionality of the map."""
dim(map::PerlinNoiseMap{N, T}) where {N, T} = return length(dims(map))

"""Compute the absolute range bound of N-dimensional Perlin noise: sqrt(N)/2."""
ran(map::PerlinNoiseMap{N, T}) where {N, T} = return sqrt(dim(map)) / 2.0

"""Compute the maximum coordinate inbound the grid for each dimension."""
function gridmax(map::PerlinNoiseMap{N, T}) where {N, T}
    # Assume each gridcell is a cube with sides of length 1.
    # Gradients correspond to gridcell corners, so there is one less
    # gradient in each dimension than there are cells.
    return size(map.gradients) .- 1
end

"""Check that the gridcell is in bounds of the noise map."""
function inbounds(map::PerlinNoiseMap{N, T}, cell::NTuple{N, Integer}) where {N, T}
    # If a dimension spans x gridpoints (gradients), then it spans x-1 gridcells.
    # Then, the valid gridcell indices are 0, 1, ..., x-2 for that dimension.
    return all(0 .<= cell .&& cell .< gridmax(map))
end

"""Find the gridcell that the point is contained in.

@post The corresponding gridcell must be inbounds of the noise map.
"""
function gridcell(map::PerlinNoiseMap{N, T}, point::NTuple{N, T}) where {N, T}
    # Assume each gridcell is a cube with sides of length 1.
    cell::NTuple{N, Integer} = Tuple( @. Integer(floor(point)) )
    # Given a gridcell with left bound lo and right bound hi.
    # Then consider a point x in bounds of the gridcell if lo <= x < hi.
    # But, for the rightmost gridcell in the grid, instead lo <= x <= hi holds.
    cell = cell .- (point .== gridmax(map))
    @assert inbounds(map, cell) "The gridcell falls outside the noise map."
    return cell
end

"""Compute the indices of the corners of the gridcell."""
function corners(cell::NTuple{N, Integer}) where {N}
    return reduce(vcat, Iterators.product([
        [cellidx, cellidx+1] for cellidx in cell
    ]...))
end

"""The smoothstep sigmoid-like smoothing function.

@pre The input must be in [0, 1].
@return A smoothed value in [0, 1].
"""
function smoothstep(x::T) where {T<:Number}
    @assert 0.0 <= x <= 1.0 "x = $x in [0, 1] does not hold."
    # Perform sigmoid-like smoothing: 3(x^2) - 2(x^3) = (3 - 2x) * x^2
    s = (3 - 2*x) * x^2
    @assert 0.0 <= s <= 1.0 "s = $s in [0, 1] does not hold."
    return s
end

"""Linearly interpolate between a and b."""
function lerp(frac::T, a::T, b::T) where {T}
    @assert 0.0 <= frac <= 1.0 "The interpolation fraction must be in [0, 1]."
    return a + frac*(b - a)
end

"""Obtain the noise value for the given point in space.

@param[in] normalize
    If true, then output a noise value in [-1, 1].
    If false, then output a noise value in [-sqrt(N)/2, sqrt(N)/2].
    Due to rounding errors, the noise may slightly fall outside the ranges.
"""
function noise(map::PerlinNoiseMap{N, T}, samplepoint::NTuple{N, T}; normalize::Bool=true) where {N,T}
    # Verify preconditions.
    @assert all(@. !isnan(samplepoint)) "The sample contains NaN values: $samplepoint."

    # Determine which grid cell the sample point falls in.
    cell::NTuple{N, Integer} = gridcell(map, samplepoint)

    # Obtain one noise value for each hypercube corner.
    # Note: Julia is 1-indexed, so increment EVERY gradient corner index.
    noise::Vector{T} = [
        dot(map.gradients[(corner .+ 1)...], samplepoint .- corner)
        for corner in corners(cell)
    ]

    # Linearly interpolate the noise value sequence.
    # The i-th iteration lerps the i-th dimension.
    for dim in 1:N
        # Each iteration halves the length of the sequence.
        @assert length(noise) == 2^(N-dim+1)

        # Determine the inerpolation fraction along the current dimension.
        # Make the linear interpolation sigmoid-like using smootstep.
        t = smoothstep(samplepoint[dim] - cell[dim])

        # Compute all interpolations along the current dimension.
        noise = [
            lerp(t, noise[i], noise[i+1])
            for i in 1:2:length(noise)
        ]
    end

    # Finally, extract the final, lerped noise value.
    @assert length(noise) == 1 "Lerping failed, too many noise values remain."
    # Assert the theoretical bounds of all generated noise values.
    @assert abs(noise[1]) <= ran(map)+1e-6 "Noise ($(noise[1])) not in valid range."
    return normalize ? noise[1] / ran(map) : noise[1]
end
