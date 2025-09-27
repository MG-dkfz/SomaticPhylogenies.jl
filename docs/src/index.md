# SomaticPhylogenies

`SomaticPhylogenies` provides a framework to simulate the accumulation and dispersion of somatic mutations across a cell population over time.
Currently, two elementary stochastic processes are implemented: [`BirthDeathProcess`](@ref) and [`MoranProcess`](@ref).
Both generate a [`Phylogeny`](@ref), which contains the full mutation history of the cell population.
Various functions exists to analyze or further process the `Phylogeny`.
For example, [`WholeGenomeSequencing`](@ref) may be simulated, [`subsampling`](@ref), or the [`MRCA`](@ref) may be determined.
For a complete list of implemented types and methods, see [`References`](@ref).

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

To run the simulation, we first have to decide how and with which rate mutations are accumulated.
The how is addressed by the abtract type [`MutationMode`](@ref); we opt for [`ReplicationCounter`](@ref).
In this case, mutations occur with every cell division, and a [`MutationRate`](@ref) generating `2` mutations per cell division is specified as follows:

```julia
μ = MutationRate(ReplicationCounter, 2)
```

Next, we have to specify the initial number of cells `N0` (defaults to `1`) and the initial time point `t0` (defaults to `0.0`).
Using the defaults, a `Phylogeny` is simulated by calling `process` with `μ`:

```julia
phylo = process(μ)
```

Next, we want to sequence the cell population.
For this, we define [`WholeGenomeSequencing`](@ref).
We have to provide the sequencing `depth` and (optionally) the bin edges for the variant allele frequencies.
Furthermore, the [`SamplingMode`](@ref) must be specified; we opt for [`PoissonBinomialSampling`](@ref).

```julia
depth = 90
vaf_edges = 0.05:0.05:1
wgs = WholeGenomeSequencing(PoissonBinomialSampling, vaf_edges, depth; read_min=3)
seq_res = wgs(phylo)
```

Finally, we can visualize the [`SequencingResult`](@ref) with [`vafhistogram`](@ref) and [`vafogram`](@ref).

```julia
using Plots

p1 = vafhistogram(seq_res)
p2 = vafogram(seq_res)
plot(p1, p2)
```