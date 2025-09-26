# References

## Cells and lineages
```@docs
Cell
exists
exists_at
CellLineage
celltime
cellgeneration
num_divisions(lineage::CellLineage)
cell_at
celllineage_at
```

## Mutations
```@docs
MutationMode
TimeCounter
ReplicationCounter
MutGenProcess
DeterministicProcess
PoissonProcess
MutationRate
mutgenprocess(μ::MutationRate{MM, MGP}) where {MM, MGP}
mutation_mode(μ::MutationRate{MM, MGP}) where {MM, MGP}
```

## Phylogeny
```@docs
Phylogeny
mutation_mode(phylo::Phylogeny{M}) where M <: MutationMode
mutgenprocess(phylo::Phylogeny)
founder_cells
population_size
subpopulation_sizes(phylo::Phylogeny)
cellgenerations
num_divisions(phylo::Phylogeny)
num_mutations(lineage::CellLineage)
lineage_vector
cell_vector(phylo::Phylogeny)
cell_count(phylo::Phylogeny)
clone_sizes
subsampling
subtree
select_lineage
site_frequency_spectrum
MRCA
leading_subclones
```

## Processes
```@docs
AbstractProcess
BirthDeathProcess
MoranProcess
CompositeProcess
```

## Whole-genome sequencing
```@docs
SamplingMode
BinomialSampling
PoissonBinomialSampling
BetaSampling
WholeGenomeSequencing
SequencingResult
sampling_mode
vaf_min(wgs::WholeGenomeSequencing)
vaf_step(wgs::WholeGenomeSequencing)
cumulative_counts
vafhistogram
vafogram
```
