
module SomaticPhylogenies

using Distributions
using StatsBase
using RecipesBase

include("Phylogeny.jl")
include("Processes.jl")
include("WholeGenomeSequencing.jl")

# Phylogeny.jl
export
    Cell,
    CellLineage,
    MutationMode,
    ReplicationCounter,
    TimeCounter,
    MutGenProcess,
    DeterministicProcess,
    PoissonProcess,
    MutationRate,
    Phylogeny;
export
    exists,
    exists_at,
    celltime,
    cellgeneration,
    cell_at,
    celllineage_at,
    mutation_mode,
    mutgenprocess,
    founder_cells,
    population_size,
    subpopulation_sizes,
    cellgenerations,
    num_divisions,
    num_mutations,
    lineage_vector,
    cell_vector,
    cell_count,
    clone_sizes,
    subsampling,
    subtree,
    select_lineage,
    site_frequency_spectrum,
    MRCA,
    leading_subclones;

# Processes.jl
export
    AbstractProcess,
    BirthDeathProcess,
    MoranProcess,
    CompositeProcess;

# WholeGenomeSequencing.jl
export
    SamplingMode,
    BetaSampling,
    BinomialSampling,
    PoissonBinomialSampling,
    WholeGenomeSequencing,
    SequencingResult;
export
    vaf_min,
    vaf_step,
    vafs,
    sampling_mode,
    cumulative_counts,
    vafhistogram,
    vafhistogram!,
    vafogram,
    vafogram!;
end
