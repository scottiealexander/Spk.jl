# Spk.jl
Julia implementation (and extension) of my older Matlab **spk** analysis package.

## Dependencies
* JSON
* Optim

Dependencies can be installed with:
```julia
import Pkg
Pkg.add(["JSON", "Optim"])
```

## Usage
The simplest way to use this package (which consists of 4 julia modules), is to modify your `startup.jl` file, which can be found in:
```julia
startuppath = joinpath(homedir(), ".julia", "config", "startup.jl")
```

Open that file and add:
```julia
spkpath = "<PATH_TO_DOWNLOADED_Spk.jl_FOLDER>"
push!(LOAD_PATH, spkpath)
```

Then (e.g.):
```julia
using SpkCore
?psth
```
