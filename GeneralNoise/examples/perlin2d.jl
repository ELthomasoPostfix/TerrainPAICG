include("../src/GeneralNoise.jl")
using .GeneralNoise

using Plots
using Random
# Prevent flaky examples by setting an RNG seed.
Random.seed!(0)

# Generate a 2D map of gradients.
map = PerlinNoiseMap{Float64}((10, 10))

# Sample noise and plot it as a heatmap.
p = image2d(map, 100, 100)
savefig(p, (@__DIR__) * "/perlin2d-plain.png")

# Sample noise and plot it as a heatmap.
# Also plot the gradients as a quiver plot.
p = image2d(map, 100, 100, pquiver=true)
savefig(p, (@__DIR__) * "/perlin2d-quiver.png")

# Sample noise and plot it as a heatmap.
# Also plot the gradients as a quiver plot.
# Also plot the sample points as a scatter plot.
p = image2d(map, 100, 100, pquiver=true, pscatter=true)
savefig(p, (@__DIR__) * "/perlin2d-quiver-scatter.png")

