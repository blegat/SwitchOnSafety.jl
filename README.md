# Switch On Safety (SOS)

| **Documentation** | **Build Status** |
|:-----------------:|:----------------:|
| [![][docs-stable-img]][docs-stable-url] | [![Build Status][build-img]][build-url] [![Build Status][winbuild-img]][winbuild-url] |
| [![][docs-latest-img]][docs-latest-url] | [![Coveralls branch][coveralls-img]][coveralls-url] [![Codecov branch][codecov-img]][codecov-url] |

This packages implements methods for computing invariant sets using [Sum Of Squares Programming](https://github.com/JuliaOpt/SumOfSquares.jl).
It supports:
* Systems defined in [MathematicalSystems.jl](https://github.com/JuliaReach/MathematicalSystems.jl).
* Hybrid Systems defined in [HybridSystems.jl](https://github.com/blegat/HybridSystems.jl).

It also includes utilities for approximation the [Joint Spectral Radius](https://link.springer.com/book/10.1007%2F978-3-540-95980-9).

## Installation

The package currently requires Julia v1.0, you can download it [here](https://julialang.org/downloads/).
It also needs the latest unreleased versions of JuMP, PolyJuMP and SumOfSquares so once Julia is installed,
launch the REPL an type
```julia
] add JuMP#master
] add PolyJuMP#master
] add SumOfSquares#master
] add https://github.com/blegat/SwitchOnSafety.jl.git
```

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-stable-url]: https://blegat.github.io/SwitchOnSafety.jl/stable
[docs-latest-url]: https://blegat.github.io/SwitchOnSafety.jl/latest

[build-img]: https://travis-ci.org/blegat/SwitchOnSafety.jl.svg?branch=master
[build-url]: https://travis-ci.org/blegat/SwitchOnSafety.jl
[winbuild-img]: https://ci.appveyor.com/api/projects/status/afm00qxhh89f8drm/branch/master?svg=true
[winbuild-url]: https://ci.appveyor.com/project/blegat/switchedsystems-jl/branch/master
[coveralls-img]: https://coveralls.io/repos/blegat/SwitchOnSafety.jl/badge.svg?branch=master&service=github
[coveralls-url]: https://coveralls.io/github/blegat/SwitchOnSafety.jl?branch=master
[codecov-img]: http://codecov.io/github/blegat/SwitchOnSafety.jl/coverage.svg?branch=master
[codecov-url]: http://codecov.io/github/blegat/SwitchOnSafety.jl?branch=master
