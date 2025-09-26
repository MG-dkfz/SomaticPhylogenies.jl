# SomaticPhylogenies

`SomaticPhylogenies` provides a framework to simulate the accumulation and dispersion of somatic mutations across a cell population over time.
Mutation accumulation relies on three types: (i) a single [`Cell`](@ref), (ii) a single [`CellLineage`](@ref), abstracted as a collection of `Cell`s, and (iii) [`Phylogeny`](@ref), a collection of lineages that result from some stochastic process ([`AbstractProcess`](@ref)).
Implemented processes are [`BirthDeathProcess`](@ref), [`MoranProcess`](@ref) or some combination thereof via [`CompositeProcess`](@ref).
`Phylogenys` also implements [`MutationMode`](@ref), which allows control over the principal process governing mutation accumulation: [`ReplicationCounter`](@ref) and [`TimeCounter`](@ref).
The package further implements [`WholeGenomeSequencing`](@ref) and [`SequencingResult`](@ref) to simulate sequencing of `Phylogeny`.

See also: [`References`](@ref)

## Installation

tbd

## Basic usage

In this example, we want to sequence a population of cells.
The first step is to specify a stochastic process.
We choose a [`BirthDeathProcess`](@ref) with division rate `λ = 1` and loss rate `δ = 0.1` that lasts for a time period of `Δt = 10`.

```julia
using Phylogeny

λ = 1    # division rate
δ = 0.1  # loss rate
process = BirthDeathProcess(λ, δ; Δt = 10)
```

Next, we run the simulation.
We have to provide the initial number of cells `N0` (defaults to `1`) and the initial time point `t0` (defaults to `0.0`).
Furthermore, the [`MutationMode`](@ref) must be specified; we opt for [`ReplicationCounter`](@ref).

```julia
N0 = 1
t0 = 0.0
mut_share = process(ReplicationCounter, N0, t0)
# equivalently:
# mut_share = process(ReplicationCounter)
```

The object `mut_share` is of type [`Phylogeny`](@ref), which holds all relevant information about mutation accumulation during the `process`.
The most important field is `mut_share.cells`, which is a `Vector` of [`Cell`](@ref)s, where each `Cell` is itself a Vector of [`Cell`](@ref)s.
`mut_share.cells` holds all information about the dispersion of `Mutation`s across the population.

Next, we want to sequence the cell population.
For this, we define [`WholeGenomeSequencing`](@ref).
We have to provide the sequencing `depth` and (optionally) the bin edges for the variant allele frequencies.
Furthermore, the [`SamplingMode`](@ref) must be specified; we opt for [`PoissonBinomialSampling`](@ref).

```julia
depth = 90
vaf_edges = 0.05:0.025:1
wgs = WholeGenomeSequencing(PoissonBinomialSampling, vaf_edges, depth; read_min=3)
```

Before we can simulate whole-genome sequencing, we have to specify how to determine the number of variants a `Mutation` encodes for.
The different ways of doing so for `ReplicationCounter`s and `TimeCounter`s are described in the documentation of [`num_mutations`](@ref).
We assume that the number of variants per cell division follows a Poisson distribution with expectation value `μ = 1.3`.

```julia
using Distributions

μ = 1.3 # average number of mutations per cell division
seq_res = wgs(mut_share, Poisson(μ))
```

Finally, we can visualize the [`SequencingResult`](@ref) with [`vafhistogram`](@ref) and [`vafogram`](@ref).

```julia
using Plots

p1 = vafhistogram(seq_res)
p2 = vafogram(seq_res)
plot(p1, p2)
```