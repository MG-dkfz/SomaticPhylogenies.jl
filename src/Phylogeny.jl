
"""
    Cell

Represents the basic element of a `CellLineage`.

## Fields
* `uid::T`        : Unique identifier of `Cell`.
* `generation::T` : Position of `Cell` within a `CellLineage`.
* `t::F`          : Time point of the `Cell`s birth.
* `mutations::M`  : Number of `Cell`'s private mutations.

## Base methods
* `==` : compares `uid`

See also: [`CellLineage`](@ref), [`MutationMode`](@ref), [`MutationRate`](@ref)
"""
mutable struct Cell{I,F}
    const uid::I
    const generation::I
    const t::F
    mutations::F
end

Cell(uid::Integer) = Cell(uid, 0.0)
Cell(uid::Integer, t::Real) = Cell(uid, t, 0.0)
Cell(uid::Integer, t::Real, mutations::Real) = Cell(uid, 1, t, mutations)
function Cell(uid::Integer, generation::Integer, t::Real, mutations::Real)
    uid < zero(uid) && error("uid must not be negative")
    Cell(promote(uid, generation)..., promote(float(t), mutations)...)
end

# struct Cell{I<:Integer, F<:AbstractFloat, M<:Real}
#     uid::I
#     t::F
#     generation::I
#     mutations::M
# end

# Cell() = Cell(0, 0.0, 0, 0.0)
# Cell(uid::Integer) = Cell(uid, 0.0)
# Cell(uid::Integer, t::Real) = Cell(uid, t, 0.0)
# Cell(uid::Integer, t::Real, mutations::Real) = Cell(uid, t, 1, mutations)
# function Cell(uid::Integer, t::Real, generation::Integer, mutations::Real)
#     uid < zero(uid) && error("uid must not be negative")
#     uid, generation = promote(uid, generation)
#     Cell(uid, float(t), generation, float(mutations))
# end


nullcell() = Cell(0, 0.0, 0, 0.0)
isnullcell(cell::Cell) = cell.uid == 0

function Base.show(io::IO, cell::Cell)
    print(io, "Cell($(cell.uid), k = $(cell.generation), t = $(cell.t))")
end
function Base.show(io::IO, ::MIME"text/plain", cell::Cell)
    s = string(typeof(cell))
    s *= "\n  * uid        : $(cell.uid)"
    s *= "\n  * generation : $(cell.generation)"
    s *= "\n  * time       : $(cell.t)"
    s *= "\n  * mutations  : $(cell.mutations)"
    print(io, s)
end

Base.:(==)(c1::Cell, c2::Cell) = c1.uid == c2.uid
Base.length(c::Cell) = 1
Base.iterate(c::Cell, k=1) = nothing
Base.broadcastable(c::Cell) = Ref(m)


"""
    exists(cell::Cell) -> Bool

Check whether `Cell` exists; returns `false` if `cell.uid == 0`.
"""
exists(cell::Cell) = cell.uid > 0

"""
    exists_at(t::Real, ::Cell) -> Bool

Check whether `Cell` exists at time `t`.
"""
exists_at(t, cell::Cell) = cell.t ≤ t

###################################################################################################

"""
    const CellLineage = Vector{Cell}

Time-continuous sequence of `Cell`s.

See also: [`Cell`](@ref)
"""
const CellLineage = Vector{<:Cell}

"""
    celltime(::CellLineage) -> Real

Return the time point of the last cell division along `CellLineage`.
"""
celltime(lineage::CellLineage) = length(lineage) > 0 ? lineage[end].t : nothing

"""
    cellgeneration(::CellLineage) -> Int

Return the length of `CellLineage`.
"""
cellgeneration(lineage::CellLineage) = length(lineage)

"""
    num_divisions(::CellLineage) -> Int

Return the number of division along `CellLineage`.
"""
num_divisions(lineage::CellLineage) = cellgenerations(lineage) - 1

"""
    cell_at(t, ::CellLineage) -> Cell

Return the `Cell` that existed at time `t` in `CellLineage`.
"""
function cell_at(t, lineage::CellLineage)
    cell = lineage[1]
    for _cell in lineage
        if exists_at(t, _cell)
            cell = _cell
        else
            break
        end
    end
    return cell
end

"""
    celllineage_at(t, lineage::CellLineage) -> CellLineage

Return the `CellLineage` that existed at time `t` in `lineage`.
"""
function celllineage_at(t, lineage::CellLineage)
    for gen in eachindex(lineage)
        if !exists_at(t, lineage[gen])
            return lineage[1:gen-1]
        end
    end
end

###################################################################################################

"""
    const Population = Vector{<:CellLineage}

See also: [`Cell`](@ref), [`CellLineage`](@ref)
"""
const Population = Vector{<:CellLineage}

initial_population(N0, t0) = [[Cell(k, t0, N0[k])] for k in eachindex(N0)]

###################################################################################################

"""
    MutationMode

Abstract type that specifies (i) how daughter cells inherit mutations from the mother cell, and (ii) how a `Mutation` translates to an actual number of mutations.

It is either [`TimeCounter`](@ref) or [`ReplicationCounter`](@ref).
"""
abstract type MutationMode end

"""
    TimeCounter <: MutationMode

Describes that the number of newly arising mutations in a cell increases linearly with the cell's age. Together with the already present mutations, these new mutations are passed on to both daughter cells.

See also: [`Cell`](@ref), [`MutationRate`](@ref)
"""
struct TimeCounter <: MutationMode end

"""
    ReplicationCounter <: MutationMode

Describes that new mutations occur with each replication event.

See also: [`MutationRate`](@ref)
"""
struct ReplicationCounter <: MutationMode end

###################################################################################################

"""
    MutGenProcess

Abstract type that specifies whether or not stochastic fluctuations in the number of newly arising mutations in a cell should be taken into account.

It is either `DeterministicProcess` or `PoissonProcess`.

See [`MutationRate`](@ref).
"""
abstract type MutGenProcess end

"""
    DeterministicProcess <: MutGenProcess

See [`MutationRate`](@ref).
"""
struct DeterministicProcess <: MutGenProcess end

"""
    PoissonProcess <: MutGenProcess

See [`MutationRate`](@ref).
"""
struct PoissonProcess <: MutGenProcess end

"""
    MutationRate{MM<:MutationMode, MGP<:MutGenProcess}

Determines mutation accumulation along a `CellLineage`.

The two parameters `MutationMode` and `MutGenProcess` specify how the single field `μ` of `MutationRate` is interpreted.

* `ReplicationCounter` : `μ` has unit **mutations per cell division** and counts for example the (average) number replication errors. Upon each cell division, each daughter cell inherits all variants from its mother plus a number of novel, private mutations. For `DeterministicProcess` this number is given by `μ`, for `PoissonProcess` `μ` serves as expectation value for a poisson distrution.
* `TimeCounter` :  `μ` has unit **mutations per cell and time** and counts for example the (average) number of DNA repair errors in between cell divisions. Since these mutations occur before the mother cell divides, both daughter cells inherit the same mutations, which comprise the mutations the mother cell inherited from its mother plus a number of novel mutations, which depend linearly on the lifespan of the mother cell, `Δt`, where `μ` serves as proportionality constant (mutation clock). For `DeterministicProcess` this number is `μ Δt`, while for `PoissonProcess` `μ Δt` serves as expectation value for a poisson distrution.

## Constructor
    MutationRate(
        [MM::Type{<:MutationMode} = ReplicationCounter,]
        μ::Real;
        poisson::Bool = false,
    )

Construct `MutationRate` with value `μ`. The default `MutationMode` is `ReplicationCounter` and `PoissonProcess` may be specified by setting `poisson = true`.

See also: [`MutationMode`](@ref), [`MutGenProcess`](@ref), [`CellLineage`](@ref)
"""
struct MutationRate{MM<:MutationMode, MGP<:MutGenProcess}
    μ::Union{Real, Poisson}
end

function MutationRate(
    μ::Real,
    MM::Type{<:MutationMode}=ReplicationCounter;
    poisson=false,
    )
    if poisson
        return MutationRate{MM, PoissonProcess}(Poisson(μ))
    else
        return MutationRate{MM, DeterministicProcess}(float(μ))
    end
end

function (μ::MutationRate{ReplicationCounter, DeterministicProcess})(cell::Cell, t)
    μ.μ
end
function (μ::MutationRate{TimeCounter, DeterministicProcess})(cell::Cell, t)
    μ.μ * (t - cell.t)
end

function (μ::MutationRate{ReplicationCounter, PoissonProcess})(cell::Cell, t)
    float(rand(μ.μ))
end
function (μ::MutationRate{TimeCounter, PoissonProcess})(cell::Cell, t)
    float(rand(Poisson(μ.μ.λ * (t - cell.t))))
end

"""
    mutgenprocess(::MutationRate) -> MutGenProcess

Return the `MutGenProcess` of `MutationRate`.
"""
mutgenprocess(μ::MutationRate{MM, MGP}) where {MM, MGP} = MGP

"""
    mutation_mode(::MutationRate) -> MutationMode

Return the `MutationMode` of `MutationRate`.
"""
mutation_mode(μ::MutationRate{MM, MGP}) where {MM, MGP} = MM

#--------------------------------------------------------------------------------------------------

"""
    Phylogeny

`Phylogeny` of a cell population at time `t`. Each member of the population is represented by a `CellLineage` holding ancestral information in terms of `Cell`s. A `Cell` defines a `uid` and contains information about the time of its birth, the number of its private mutations, and its position within a `CellLineage`. Population members descending from the same cell all have that `Cell` at the same position in their `CellLineage`.

The `Phylogeny`-generating dynamics are simulated with `AbstractProcess`es  .

## Fields
* `lineages::Vector{CellLineage}`
* `t::Real` : The time point of observation.
* `μ::MutationRate`
* `uid_counter::Int` : The global uid counter, which ensures that new mutations in a `Cell` are uniquely identified (infinte-site hyothesis).

See also: [`Cell`](@ref), [`CellLineage`](@ref), [`MutationRate`](@ref), [`AbstractProcess`](@ref)

---

## Constructors
    Phylogeny(μ::MutationRate[, N0::Integer, t0])

Create an initial `Phylogeny` with `N0` (default `1`) somatically unmutated founder cells with initial time `t0` (default `0.0`).

    Phylogeny(μ::MutationRate, N0::Vector{<:Real}[, t0])

In order to equip the founder cells with somatic mutations, a vector `N0` may be passed where `length(N0)` is the number of founder cells and the elements are the number of founder mutations.

See also: [`MutationRate`](@ref)
"""
struct Phylogeny{K,L,M}
    lineages::NTuple{K, L}
    t::AbstractFloat
    t0::AbstractFloat
    mutation_rate::NTuple{K, M}
    uid_counter::Integer
end

# Constructors
# -------------

# i) initial Phylogeny{1} from MutationRate without founder mutations
function Phylogeny(
    μ::MutationRate,
    N0::Integer=1;
    t0::Real=0.0
    )
    Phylogeny(μ, zeros(N0); t0=t0)
end

# ii) initial Phylogeny{1} from MutationRate with founder mutations
function Phylogeny(
    μ::MutationRate,
    N0::Vector{<:Real};
    t0::Real=0.0)
    
    Phylogeny(initial_population(N0, t0), t0, t0, μ, length(N0))
end

# iii) Phylogeny{1} from population and MutationRate
function Phylogeny(
    lineages::Population,
    t0, t,
    μ::MutationRate,
    uid_counter
    )
    
    Phylogeny((lineages, ), t0, t, (μ, ), uid_counter)
end

# iv) initial Phylogeny{K} from MutationRate without founder mutations
function Phylogeny(
    μ::Vector{<:MutationRate},
    N0::Integer=1;
    t0::Real=0.0
    )
    Phylogeny(μ, zeros(N0); t0=t0)
end

# v) initial Phylogeny{K} from MutationRate with founder mutations
function Phylogeny(
    μ::Vector{<:MutationRate},
    N0::Vector{<:Real};
    t0::Real=0.0)
    
    Phylogeny(initial_population(N0, t0), t0, t0, μ, length(N0))
end

# vi) Phylogeny{K} from population and K MutationRates
function Phylogeny(
    pop::Population,
    t0, t,
    μ::Vector{<:MutationRate},
    uid_counter
    )
    populations = tuple(pop, fill(empty_population(), length(μ)-1)...)
    Phylogeny(populations, t0, t, tuple(μ), uid_counter)
end

# iv) updated Phylogeny{1}
function Phylogeny(phylo::Phylogeny, lineages, t, uid_counter)
    Phylogeny(lineages, promote(phylo.t0, t)..., phylo.mutation_rate, uid_counter)
end


function (phylo::Phylogeny{ReplicationCounter})(lineage, uid, t)
    mother = lineage[end]
    daughter1 = push!(
        lineage,
        Cell(uid+1, t, mother.generation+1, phylo.mutation_rate(mother, t))
        )
    daughter2 = copy(daughter1)
    daughter2[end] = Cell(uid+2, t, mother.generation+1, phylo.mutation_rate(mother, t))
    return daughter1, daughter2, uid + 2
end

function (phylo::Phylogeny{TimeCounter})(lineage, uid, t)
    mother = lineage[end]
    mother.mutations = phylo.mutation_rate(mother, t)
    daughter1 = push!(
        lineage,
        Cell(uid+1, t, mother.generation+1, 0.0)
        )
    # daughter1 = push!(
    #     lineage,
    #     Cell(uid+1, t, mother.generation+1, phylo.mutation_rate(mother, t))
    #     )
    return daughter1, copy(daughter1), uid+1
end

#--------------------------------------------------------------------------------------------------

function Base.show(io::IO, phylo::Phylogeny{1})
    s  = "$(typeof(phylo))"
    s *= "\n  * Time             : $(phylo.t)"
    s *= "\n  * Population size  : $(length(phylo.lineages))"
    s *= "\n  * Mutation counter : $(phylo.uid_counter)"
    print(io, s)
end

Base.length(phylo::Phylogeny{1}) = length(phylo.lineages[1])
Base.firstindex(phylo::Phylogeny{1}) = 1
Base.lastindex(phylo::Phylogeny{1}) = length(phylo)
function Base.getindex(phylo::Phylogeny{1}, i::Int)
    1 ≤ i ≤ length(phylo) || throw(BoundsError(phylo, i))
    return phylo.lineages[1][i]
end
Base.getindex(phylo::Phylogeny{1}, i::Number) = phylo[convert(Int, i)]
Base.getindex(phylo::Phylogeny{1}, I) = [phylo[i] for i in I]
Base.iterate(phylo::Phylogeny{1}, k=1) = k > length(phylo) ? nothing : (phylo[k], k+1)
Base.broadcastable(phylo::Phylogeny) = Ref(phylo)


###################################################################################################

"""
    mutation_mode(::Phylogeny) -> MutationMode

Return the `MutationMode` of `Phylogeny`.
"""
mutation_mode(phylo::Phylogeny) = mutation_mode(phylo.mutation_rate)

"""
    mutgenprocess(::Phylogeny) -> MutGenProcess

Return the `MutGenProcess` of `Phylogeny.mutation_rate`.
"""
mutgenprocess(phylo::Phylogeny) = mutgenprocess(phylo.mutation_rate)

"""
    founder_cells(::Phylogeny) -> Vector{Cell}

Return the founder cells of `Phylogeny`.
"""
function founder_cells(phylo::Phylogeny)
    unique(cell_at(phylo.t0, lin) for lin in phylo)
end

"""
    population_size(::Phylogeny) -> Int

Return the population size of `Phylogeny`.
"""
population_size(phylo::Phylogeny) = length(phylo)

"""
    subpopulation_sizes(phylo::Phylogeny) -> Dict{Cell, Int}

Return the number of cells that derive from distinct founder cells.

    subpopulation_sizes(generation::Integer, phylo::Phylogeny) -> Dict{Cell, Int}

Instead of founder cells, a `generation` may be provided whose `Cell`s define the sub-populations.

    subpopulation_sizes(t::AbstractFloat, phylo::Phylogeny) -> Dict{Cell, Int}

Instead of `generation`, a time point `t` may be provided where `Cell`s that existed at `t` define the sub-populations.
"""
function subpopulation_sizes(phylo::Phylogeny)
    subpopulation_sizes(1, phylo)
end
function subpopulation_sizes(generation::Integer, phylo::Phylogeny)
    if generation > minimum(cellgenerations(phylo))
        throw(DomainError(generation, "not all `lineages` have reached generation $(generation)"))
    end
    # count `Mutation`s
    countmap([lineage[generation] for lineage in phylo])
end
function subpopulation_sizes(t::AbstractFloat, phylo::Phylogeny)
    if t > phylo.t
        throw(DomainError(t, "time `t` greater than time of observation"))
    end
    # count `Mutation`s
    countmap([cell_at(t, lineage) for lineage in phylo])
end

"""
    cellgenerations(phylo::Phylogeny) -> Vector{Int}

Return the lineage lengths of `Phylogeny`.
"""
cellgenerations(phylo::Phylogeny) = cellgeneration.(phylo.lineages)

"""
    num_divisions(phylo::Phylogeny) -> Vector{Int}

Return the number of divisions each lineage of `Phylogeny` underwent.
"""
num_divisions(phylo::Phylogeny) = cellgenerations(phylo) .- 1

"""
    num_mutations(::CellLineage) -> Real

Return the number of mutations acquired along `CellLineage`.

    num_mutations(::Phylogeny) -> Vector{<:Real}

Return the number of mutations acquired along each `CellLineage` of `Phylogeny`.
"""
function num_mutations(lin::CellLineage)
    sum(cell.mutations for cell in lin)
end
function num_mutations(phylo::Phylogeny)
    [num_mutations(lin) for lin in phylo]
end

"""
    lineage_vector(::Phylogeny) -> Vector{Cell}

Return concatanation of all `CellLineage`s in `Phylogeny`.
"""
lineage_vector(phylo::Phylogeny) = vcat(phylo.lineages...)

"""
    cell_vector(::Phylogeny) -> Vector{Cell}
    
Return all `Cell`s that ever occurred in `Phylogeny`.

    cell_vector(generation::Integer, ::Phylogeny) -> Vector{Cell}

Return all `Cell`s of `generation` or greater that ever occurred in `Phylogeny`.
"""
cell_vector(phylo::Phylogeny) = unique(lineage_vector(phylo))
cell_vector(generation::Integer, phylo::Phylogeny) = begin
    unique([lineage[generation] for lineage in phylo if length(lineage) ≥ generation])
end

"""
    cell_count(::Phylogeny [, n_min]) -> Int

Return the total number of distinct `Cell`s that ever occurred in Phylogeny.

If `n_min` is specified, only `Cell`s that occur in at least `n_min` `CellLineage`s are considered.
"""
cell_count(phylo::Phylogeny) = length(cell_vector(phylo))
function cell_count(phylo::Phylogeny, n_min)
    _clone_sizes = values(clone_sizes(phylo))
    return sum(_clone_sizes .>= n_min)
end

"""
    clone_sizes(::Phylogeny) -> Dict{Cell, Int}

Return the number of lineages each `Cell` occurs in.
"""
clone_sizes(phylo::Phylogeny) = countmap(lineage_vector(phylo))

"""
    subsampling(::Phylogeny, n::Integer) -> Phylogeny

Return sub-`Phylogeny` comprising `n` randomly selected lineages (without replacement).
"""
function subsampling(phylo::Phylogeny, n::Integer)
    n > population_size(phylo) && error("`n` must not exceed population size.")
    Phylogeny(
        sample(phylo.lineages, n; replace = false),
        phylo.t,
        phylo.t0,
        phylo.mutation_rate,
        phylo.uid_counter)
end

"""
    subtree(::Phylogeny, cell::Cell) -> Phylogeny

Return sub-`Phylogeny` comprising all lineages that derive from `cell`.
"""
function subtree(phylo::Phylogeny, cell::Cell)
    Phylogeny(
        [lin for lin in phylo if (cellgeneration(lin) ≥ cell.generation) && (lin[cell.generation] == cell)],
        phylo.t,
        cell.t,
        phylo.mutation_rate,
        phylo.uid_counter)
end

"""
    select_lineage(::Phylogeny) -> CellLineage

Return a randomly selected `CellLineage` from `Phylogeny`.
"""
function select_lineage(phylo::Phylogeny)
    phylo[rand(1:population_size(phylo))]
end

###################################################################################################

"""
    site_frequency_spectrum(phylo::Phylogeny) -> Vector{<:Real}

Return the site frequency spectrum, which has the same length as `population_size(phylo)`.
"""
function site_frequency_spectrum(phylo::Phylogeny)
    site_frequency_spectrum!(zeros(population_size(phylo)), phylo)
end
function site_frequency_spectrum!(sfs::Vector{<:Real}, phylo::Phylogeny)
    mut_sizes = clone_sizes(phylo)
    for (cell, clone_size) in mut_sizes
        if clone_size > 0
            sfs[clone_size] += cell.mutations
        end
    end
    return sfs
end
function site_frequency_spectrum(phylos::Vector{<:Phylogeny})
    sfs = zeros(maximum(population_size.(phylos)))
    for phylo in phylos
        site_frequency_spectrum!(sfs, phylo)
    end
    return sfs
end

###################################################################################################

"""
    MRCA(::Phylogeny) -> CellLineage

Return the `CellLineage` of the `M`ost `R`ecent `C`ommon `A`ncestor of `Phylogeny`. 

If no MRCA exists, an empty `CellLineage` is returned.
"""
function MRCA(phylo::Phylogeny)
    mrca = nothing
    
    # if first generation is not clonal -> no clonal mutations at all
    L = 1
    if !__isclonal(L, phylo)
        return CellLineage[]
    end
    
    # minimal requirement for clonal mutation:
    #   the `generation` index of a clonal `Mutation` cannot be
    #   greater than the shortest lineage in `phylo.lineages`
    R = minimum(cellgenerations(phylo))
    
    # binary search
    m = 1
    while R >= L
        m = L + floor(Int, (R - L) / 2)
        if __isclonal(m, phylo)
            L = m + 1
        else
            R = m - 1
        end
    end
    
    if __isclonal(m, phylo)
        return __clonal_mutations(m, phylo)
    else
        return __clonal_mutations(m - 1, phylo)
    end
end

function __clonal_mutations(generation, phylo)
    phylo[1][1:generation]
end

function __isclonal(generation, phylo)
    mutation = phylo[1][generation]
    candidates = [lineage[generation] for lineage in phylo]
    for candidate in candidates
        if candidate != mutation
            return false
        end
    end
    return true
end

###################################################################################################

"""
    leading_subclones(::Phylogeny, depth::Integer)

Determine the leading subclones of `Phylogeny` up to `depth`.

## Returns
`clone_structure::Vector{Vector{Pair{<:Cell, NamedTuple{(:clonesize, :Δnmut), Tuple{Int64, Float64}}}}}` : `Vector` of length `depth` where each element represents the subclones of a level, where `level = 1:depth`. Each level contains `2^(level - 1)` clones, which are represented by a `Vector` of that length. Each clone is represented by a `Pair`, where `first` is the last `Cell` along the sublineage defining the clone, and `second` is a `NamedTuple` containing the `clonesize` and the number of mutations `Δnmut` accumulated along the sublineage.

The relationship between clones across levels is as follows: The first two clones on level `L` descent from the first clone on level `L-1`, the third and fourth clones on level `L` descent from the second clone on level `L-1`, and so on. Furthermore, it is possible that a clone on some level `L` represents a single `Cell` (`clonesize = 1`), which is therefore still alive at the time of observation. In order to keep above relationship logic valid, 
"""
function leading_subclones(phylo, depth)
    
    if isempty(MRCA(phylo))
        return nothing
    end
    
    function mrca_and_subphylos(phylo, gen)
        if population_size(phylo) == 0
            return Cell() => (clonesize = 0, Δnmut = 0.0), [phylo, phylo], gen
        end
        if population_size(phylo) == 1
            cell = phylo[1][end]
            phylo = Phylogeny(phylo.mutation_rate, 0)
            return cell => (clonesize = 1, Δnmut = cell.mutations), [phylo, phylo], gen
        end
        mrca = MRCA(phylo)
        daughters = unique(lin[mrca[end].generation + 1] for lin in phylo)
        local subphylos = [subtree(phylo, daughter) for daughter in daughters]
        return (
            mrca[end] => (clonesize = population_size(phylo), Δnmut = num_mutations(mrca[gen:end])),
            subphylos, mrca[end].generation+1)
    end

    indices(level) = 2^(level-1):2^(level)-1
    mrca, _subphylos, gen = mrca_and_subphylos(phylo, 1)
    subclones = [mrca]
    subphylos = [_subphylos]
    gens = [[gen, gen]]
    for level = 2:depth
        for k = indices(level-1)
            for (subphylo, gen) in zip(subphylos[k], gens[k])
                mrca, _subphylos, gen = mrca_and_subphylos(subphylo, gen)
                push!(subclones, mrca)
                push!(subphylos, _subphylos)
                push!(gens, [gen , gen])
            end
        end
    end
    return [subclones[indices(level)] for level = 1:depth]
end

###################################################################################################