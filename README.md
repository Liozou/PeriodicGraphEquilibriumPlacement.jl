# PeriodicGraphEquilibriumPlacement

[![Build Status](https://github.com/Liozou/PeriodicGraphEquilibriumPlacement.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/Liozou/PeriodicGraphEquilibriumPlacement.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/Liozou/PeriodicGraphEquilibriumPlacement.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/Liozou/PeriodicGraphEquilibriumPlacement.jl)

A Julia package for computing the *equilibrium*, or *barycentric*, placement of vertices
of a periodic graph, as defined by [Olaf Delgado-Friedrichs and Michael O'Keeffe](https://doi.org/10.1107/S0108767303012017).
It is accessible through the `equilibrium` exported function, which returns a matrix of
rational coordinates that can be fed to the [`PeriodicGraphEmbedding`](https://liozou.github.io/PeriodicGraphEmbeddings.jl/dev/types/#PeriodicGraphEmbeddings.PeriodicGraphEmbedding-Union{Tuple{T},%20Tuple{D},%20Tuple{PeriodicGraph{D},%20AbstractMatrix{T},%20Cell}}%20where%20{D,%20T})
or
[`SortedPeriodicGraphEmbedding`](https://liozou.github.io/PeriodicGraphEmbeddings.jl/dev/types/#PeriodicGraphEmbeddings.SortedPeriodicGraphEmbedding-Union{Tuple{T},%20Tuple{D},%20Tuple{PeriodicGraph{D},%20AbstractMatrix{T}%20where%20T,%20Cell}}%20where%20{D,%20T})
methods from [PeriodicGraphEmbeddings.jl](https://github.com/Liozou/PeriodicGraphEmbeddings.jl):

```julia
julia> tbo = PeriodicGraph3D("3 1 2 0 0 0 1 3 0 0 0 1 4 0 0 0 2 5 0 0 0 2 6 0 0 0 2 7 0 0 0 3 6 0 0 1 3 8 0 0 0 3 9 0 0 0 4 6 1 0 0 4 10 0 0 0 4 11 0 0 0 5 12 0 0 0 5 13 0 0 0 7 12 1 1 -1 7 13 0 1 0 8 12 0 0 0 8 14 0 0 0 9 12 1 1 0 9 14 0 1 0 10 13 0 0 0 10 14 0 0 0 11 13 1 1 0 11 14 1 1 -1");

julia> equilibrium(tbo)
3×14 Matrix{Rational{Int64}}:
 0//1  -1//6  -1//6   1//3  -1//3  -1//3   0//1  -1//3  0//1   0//1   2//3  -2//3  -1//6  -1//6
 0//1   0//1   0//1   0//1  -1//3   0//1   1//3  -1//3  1//3  -1//3   1//3  -1//2  -1//2  -1//2
 0//1  -1//6   1//3  -1//6   0//1  -1//3  -1//3   1//3  1//3   0//1  -1//3   1//3  -1//6   1//3
```

The implementation is optimized through a custom solver specialized for the exact
resolution of sparse integer linear system through [Dixon's algorithm](https://doi.org/10.1007/bf01459082).
The solver is directly accessible through the `dixon_solve` function:

```julia
julia> A = sparse([-3 0 2 0; 0 -5 2 3; 2 2 -2 0; 0 3 0 -3]);

julia> Y = [1 1; 0 2; 1 -1; 0 0];

julia> A * dixon_solve(Val(2), A, Y) == Y
true
```

The first argument of `dixon_solve` must be `Val(size(Y)[2])` and the second must be square.

The package also exposes a `rational_solve` function which solves the same systems through
a simpler LU decomposition approach. It serves as fallback to `dixon_solve` when Dixon's
algorithm fails, but can also be used as-is with the same API. Its performance is in
general lower than `dixon_solve`, often significantly so.

See also:

- [PeriodicGraphs.jl](https://github.com/Liozou/PeriodicGraphs.jl) for the
  underlying library and the API of the `PeriodicGraph` type.
- [CrystalNets.jl](https://github.com/coudertlab/CrystalNets.jl) for a dependent package
  specialized on crystal nets.
