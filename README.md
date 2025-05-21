# TerrainPAICG
Harnessing basic Procedural Content Generation (PCG) techniques to generate image datasets as the basis to train an image-based AI model.

## Project Structure

This project blends different technologies for the sake of it, for the joy unique to creating an overengineered mess.

Thus, the project structure is as follows.
1. A C library of basic PCG techniques and utilities for image generation.
2. A Julia module that trains image-based AI models from PCG-generated datasets, by invoking the C library.

## Testing a new Package

When adding a new local module to the project, several files must be updated for the CI pipeline to incorporate it.

First, add the package to the root `/Project.toml` file.

```toml
# Project.toml
[deps]
NewPackage = "UUID"

[sources]
NewPackage = {path = "./path/to/NewPackage"}
```

Next, add the package name to the topmost `/test/runtests.jl` file, so that its tests are also run when the CI pipeline is invoked.

```julia
Pkg.test([
    ...,
    "NewPackage"
])
```

Finally, the package can be imported inside the topmost `/src/TerrainPAICG.jl` module file, so that the source code may invoke the custom module.

```julia
module TerrainPAICG
    # Import the new module; it can now be used in the code.
    using NewPackage
end
```


## Manual Testing

To manually execute the tests associated with a Julia module, use Julia's REPL.

```sh
$ julia
julia> ]
pkg> activate .
(TerrainPAICG) pkg> test
```
