
# Sampling mode ###################################################################################

"""
    SamplingMode

Abstract type that specifies how whole-genome sequencing is simulated.

It is either [`BinomialSampling`](@ref), [`PoissonBinomialSampling`](@ref), or [`BetaSampling`](@ref).
"""
abstract type SamplingMode end

"""
    BinomialSampling <: SamplingMode

The simulated read counts for variants are assumed to be binomially distributed while the coverage for each genomic site is the same (= sequencing depth).
"""
struct BinomialSampling <: SamplingMode end

"""
    PoissonBinomialSampling <: SamplingMode

The simulated read counts for variants are assumed to be binomially distributed and the coverage of a genomic site is assumed to be poisson-distributed (sequencing depth = mean(poisson)).
"""
struct PoissonBinomialSampling <: SamplingMode end

"""
    BetaSampling <: SamplingMode

Instead of simulating read counts and coverages, the variant allele frequencies are directly determined with a Beta distribution.

* As the number of variants at true allele frequency `x` becomes larger, while each individual variant follows `BinomialSampling`, the measured variant allele frequencies converge against a Beta distribution characterized by `α = x * (seq_depth - 1), β = (1 - x) * (seq_depth - 1)`.
"""
struct BetaSampling <: SamplingMode end


# Whole-genome sequencing #########################################################################

"""
    WholeGenomeSequencing{SM<:SamplingMode}

Object to simulate whole-genome sequencing.

## Fields
* `depth::Int` : The (average) coverage of genomic sites.
* `vaf_edges::Vector{Float64}` : The bin edges for variant allele frequencies
* `read_min::Int` : Only variants with `number of reads ≥ read_min` are considered in sequencing output.
---

## Constructors
    WholeGenomeSequencing(
        ::Type{<:SamplingMode},
        [vaf_edges::AbstractVector{<:AbstractFloat},]
        depth::Integer;
        read_min::Integer=1,
        read_step::Integer=1)

Create `WholeGenomeSequencing` for `SamplingMode` and (average) coverage `depth`. If `vaf_edges` is not provided, it will be specified according to `unique(push!(collect(read_min/depth : read_step/depth : 1), 1))`.

`SamplingMode` is either `BetaSampling`, `BinomialSampling`, or `PoissonBinomialSampling`.

#### Optional keyword arguments
* `read_min::Integer=1` : The minimal read count to use in sequencing output.
* `read_step::Integer=1` : The number of reads between two adjacent `vaf_edges` (only relevant if `vaf_edges` is not supplied).

See also: [`SamplingMode`](@ref)
---

## Functor
    (wgs::WholeGenomeSequencing{<:SamplingMode})(::Phylogeny) -> SequencingResult
    (wgs::WholeGenomeSequencing{<:SamplingMode})(::Vector{<:Real}) -> SequencingResult

Simulate `WholeGenomeSequencing` and return `SequencingResult`. The input is either a `Phylogeny` or an already computed site frequency spectrum.

See also: [`SamplingMode`](@ref), [`site_frequency_spectrum`](@ref)

## Returns
A `SequencingResult` comprising two fields:
* `wgs::WholeGenomeSequencing`    : The `WholeGenomeSequencing` object used for simulations.
* `vaf_histogram::Vector{<:Real}` : The simulated distribution of variant allele frequencies across `wgs.vaf_edges`.

See also: [`SequencingResult`](@ref)
"""
struct WholeGenomeSequencing{SM<:SamplingMode}
    vaf_edges
    depth
    read_min
end

# Constructors ------------------------------------------------------------------------------------

WholeGenomeSequencing(SM::Type{<:SamplingMode}, vaf_edges, depth; read_min=1, read_step=1) = begin
    WholeGenomeSequencing{SM}(
        vaf_edges,
        depth,
        read_min,
        )
end

WholeGenomeSequencing(SM::Type{<:SamplingMode}, depth; read_min=1, read_step=1) = begin
    WholeGenomeSequencing{SM}(
        unique(push!(collect(read_min / depth : read_step / depth : 1), 1)),
        depth,
        read_min,
        )
end

Base.show(io::IO, wgs::WholeGenomeSequencing) = print(io, "$(typeof(wgs))($(wgs.depth)X)")

# Functors ----------------------------------------------------------------------------------------

function __cumulative_to_hist(counts)
    counts[1:end] .- push!(counts[2:end], zero(eltype(counts)))
end

function __to_int(mutation_counts)
    n = trunc(mutation_counts)
    return Int(n + rand(Bernoulli(mutation_counts - n)))
end


function (wgs::WholeGenomeSequencing{BetaSampling})(SFS::AbstractVector)
    observed = zeros(Float64, length(wgs.vaf_edges))
    N = length(SFS)
    for (n, counts) in enumerate(SFS)
        μ = n / 2N
        α = μ * (wgs.depth - 1)
        β = (1 - μ) * (wgs.depth - 1)
        observed .+= counts .* (1 .- cdf.(Beta(α, β), wgs.vaf_edges))
    end
    return SequencingResult(wgs, __cumulative_to_hist(observed))
end

function (wgs::WholeGenomeSequencing{BinomialSampling})(SFS::AbstractVector)
    observed = zeros(Int, length(wgs.vaf_edges))
    N = length(SFS)
    for (n, counts) in enumerate(SFS)
        counts = __to_int(counts)
        reads = rand(Binomial(wgs.depth, n/2N), counts)
        vafs = reads[reads .≥ wgs.read_min] ./ wgs.depth
        observed .+= [sum(vafs .≥ vaf) for vaf in wgs.vaf_edges]
    end
    return SequencingResult(wgs, __cumulative_to_hist(observed))
end

function (wgs::WholeGenomeSequencing{PoissonBinomialSampling})(SFS::AbstractVector)
    observed = zeros(Int, length(wgs.vaf_edges))
    poi = Poisson(wgs.depth)
    N = length(SFS)
    for (n, counts) in enumerate(SFS)
        counts = __to_int(counts)
        coverages = rand(poi, counts)
        reads = [rand(Binomial(coverage, n/2N)) for coverage in coverages]
        vafs = reads[reads .≥ wgs.read_min] ./ wgs.depth
        observed .+= [sum(vafs .≥ vaf) for vaf in wgs.vaf_edges]
    end
    return SequencingResult(wgs, __cumulative_to_hist(observed))
end

function (wgs::WholeGenomeSequencing{SM})(phylo::Phylogeny) where SM <: SamplingMode
    wgs(site_frequency_spectrum(phylo))
end

# Sequencing result ###############################################################################

"""
    SequencingResult

Result of whole-genome sequencing.

### Fields
* `wgs::WholeGenomeSequencing`    : The `WholeGenomeSequencing` object used for simulations.
* `vaf_histogram::Vector{<:Real}` : The simulated number of variant allele frequencies (`vaf`s) across `wgs.vaf_edges`. Specifically, the `k`-th value counts `vaf_edges[k] ≤ vaf < vaf_edges[k+1]`, where 1 ≤ k < length(`vaf_edges`); the number of variants with `vaf = 1` is `vaf_histogram[end]`.

See also: [`WholeGenomeSequencing`](@ref)
"""
struct SequencingResult
    wgs::WholeGenomeSequencing
    vaf_histogram
end

Base.show(io::IO, sr::SequencingResult) = print(io, "SequencingResult for $(sr.wgs)")

# Plot recipes ------------------------------------------------------------------------------------

"""
    vafhistogram(::SequencingResult)
    vafhistogram!(::SequencingResult)

Plot the simulated vaf histogram.
"""
vafhistogram

"""
    vafogram(::SequencingResult)
    vafogram!(::SequencingResult)

Plot the number of variants greater-or-equal than `vaf` versus inverted `vaf`.
"""
vafogram

@userplot VafHistogram
@recipe function plot_vafhistogram(p::VafHistogram)
    sr = p.args[1]
    seriestype := :bar
    primary --> false
    bar_width --> vaf_step(sr)
    xticks --> [0, 0.25, 0.5, 0.75, 1]
    xguide --> "VAF"
    yguide --> "Number of mutations"
    x = sr.wgs.vaf_edges[1:end-1]
    y = sr.vaf_histogram[1:end-1]
    y[end] += sr.vaf_histogram[end]
    return x, y
end


@userplot Vafogram
@recipe function plot_vafogram(p::Vafogram)
    sr = p.args[1]
    seriestype := :path
    primary --> false
    xguide --> "VAF"
    xticks --> ([1, 4, 10, 20], ["1", "0.25", "0.1", "0.05"])
    yguide --> "Cumulative number\nof mutations"
    return 1 ./ sr.wgs.vaf_edges, cumulative_counts(sr)
end

# Functions #######################################################################################

"""
    sampling_mode(::WholeGenomeSequencing) -> Type{<:SamplingMode}

Return the `SamplingMode` of `WholeGenomeSequencing` instance.
"""
sampling_mode(wgs::WholeGenomeSequencing{T}) where T = T

"""
    vaf_min(::WholeGenomeSequencing) -> Float64
    vaf_min(::SequencingResult) -> Float64

Return the lowest value of `vaf_edges`.
"""
vaf_min(wgs::WholeGenomeSequencing) = wgs.vaf_edges[1]
vaf_min(sr::SequencingResult) = vaf_min(sr.wgs)

"""
    vaf_step(::WholeGenomeSequencing) -> Float64
    vaf_step(::SequencingResult) -> Float64

Return the width of the first bin of `vaf_edges`.
"""
vaf_step(wgs::WholeGenomeSequencing) = wgs.vaf_edges[2] - wgs.vaf_edges[1]
vaf_step(sr::SequencingResult) = vaf_step(sr.wgs)

"""
    cumulative_counts(::SequencingResult) -> Vector{Float64}

Return the cumulative mutation counts of `SequencingResult`.
"""
function cumulative_counts(sr::SequencingResult)
    (sum(sr.vaf_histogram) .- 
    append!([zero(eltype(sr.vaf_histogram))], cumsum(sr.vaf_histogram[1:end-1]))
    )
end

###################################################################################################
