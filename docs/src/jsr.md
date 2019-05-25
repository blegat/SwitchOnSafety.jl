# Joint Spectral Radius (JSR)

The Joint Spectral Radius (JSR) of a discrete-time switched system is the
minimal value of `γ` such that the system scaled by `γ` is asymptotically
stable.
```@docs
ScaledHybridSystem
```
The JSR is NP-hard to compute but several methods exist to
approximate this quantity.

## Gripenberg algorithm, a Branch-and-Bound approach

The following approachs searches through all the different switching sequences,
pruning some using the running estimate of the lower bound.

```@docs
gripenberg
```

## Sum-of-Squares approach

The following method computes upper bounds using
[Sum Of Squares Programming](https://github.com/JuliaOpt/SumOfSquares.jl)
```@docs
soslyap
soslyapbs
```

The infeasibility certificates computed in the binary search carried out by
[`soslyapbs`](@ref) can be used to produce cycles of high growth rate using
[`findsmp`](@ref) on the sequence produced by [`sosbuildsequence`](@ref).

```@docs
sosbuildsequence
findsmp
```

Alternatively, [`sosextractcycle`](@ref) can be used to find cycles of high
growth rate.

```@docs
sosextractcycle
```

## Polytopic approach

The following method can verify numerically that a cycle is an s.m.p. by
computing invariant polytopes. When it is not an s.m.p., it find a candidate
of higher growth rate.

```@docs`
invariant_polytopes
```
