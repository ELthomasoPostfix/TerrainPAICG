include("TerrainPAICG.jl")
using .TerrainPAICG

# Prevent flaky tests by setting an RNG seed.
using Random
Random.seed!(0)

map = PerlinNoiseMap{2,Float64}((2, 2))
println(typeof(map.gradients)); println(map.gradients); println()
println("noise: $(noise(map, (0.5,0.5)))")



# # TODO: If a 2D noise map represents a terrain topology, then
# #          a 3D noise map represents the evolution of a topoly
# #       throughout time. I.e. if you take the third
# #       noise map dimension to represent the time axis, and assign
# #       it a suitable small time step, then you can animate the
# #       evolution of topoly through time (make a gif of this).


p = image2d(map, 400, 400, pquiver=true)

display(p)
println("Press ENTER to continue.")
readline()
