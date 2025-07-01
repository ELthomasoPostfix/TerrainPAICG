"""A library containing generalized coherent noise algorithms."""
module GeneralNoise

    using Plots
    using LinearAlgebra

    include("perlin.jl")
    include("plot.jl")

    export PerlinNoiseMap, noise
    export image1d, image2d
end