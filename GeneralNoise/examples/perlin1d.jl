include("../src/GeneralNoise.jl")
using .GeneralNoise

using Plots
using Random
# Prevent flaky examples by setting an RNG seed.
Random.seed!(0)

# Generate a 1D map of gradients.
map = PerlinNoiseMap{Float64}((20,))

# Sample noise and plot it as a curve.
p = image1d(map, 500)
savefig(p, (@__DIR__) * "/perlin1d-plain.png")

# Sample noise and plot it as a curve.
# Also plot the gradients as a quiver plot.
p = image1d(map, 500, pquiver=true)
savefig(p, (@__DIR__) * "/perlin1d-quiver.png")

# Sample noise and plot it as a curve.
# Also plot the gradients as a quiver plot.
# Also plot the sample points as a scatter plot.
p = image1d(map, 500, pquiver=true, pscatter=true)
savefig(p, (@__DIR__) * "/perlin1d-quiver-scatter.png")

