
"""
    AbstractProcess

Abstract supertype of all population-dynamics processes.

There are three commonalties for the concrete types
1. A functor is implemented
2. The functor accepts a `Phylogeny` instance as input
3. The functor returns a `Phylogeny` instance
"""
abstract type AbstractProcess end

###################################################################################################

"""
    BirthDeathProcess{F,I} <: AbstractProcess

`BirthDeathProcess` with division rate `λ` and loss rate `δ`.

`BirthDeathProcess` has three termination conditions: duration `Δt`, minimal population size `N_min`, and maximal population size `N_max`.

---

## Constructors
    BirthDeathProcess(
        λ::Real,
        δ::Real;
        Δt::Real = Inf,
        N_min::Real = 0,
        N_max::Real = Inf)
    BirthDeathProcess(
        λ::Real = 1.0;
        Δt::Real = Inf,
        N_min::Real = 0,
        N_max::Real = Inf)

Construct `BirthDeathProcess`. If `δ` is not provided, it is set to `0.0` (pure birth process).
---

## Functors
    (::BirthDeathProcess)(μ::MutationRate [, N0, t0::Real]) -> Phylogeny
    (::BirthDeathProcess)(phylo::Phylogeny) -> Phylogeny

Simulate `BirthDeathProcess` and return `Phylogeny`.
"""
struct BirthDeathProcess{F,I} <: AbstractProcess
    λ::F
    δ::F
    Δt::F
    N_min::I
    N_max::I
end

function BirthDeathProcess(λ, δ; Δt=Inf, N_min=0, N_max=Inf)
    λ, δ, Δt, _ = promote(λ, δ, Δt, 1.0)
    N_min, N_max = promote(N_min, N_max)
    BirthDeathProcess(λ, δ, Δt, N_min, N_max)
end

BirthDeathProcess(λ=1.0; kwargs...) = BirthDeathProcess(λ, 0.0; kwargs...)

function Base.show(io::IO, process::BirthDeathProcess)
    s = "$(typeof(process))"
    s *= "\n  * Division rate λ   : $(process.λ)"
    s *= "\n  * Loss rate δ       : $(process.δ)"
    s *= "\n  * Duration Δt       : $(process.Δt)"
    s *= "\n  * Lower bound N_min : $(process.N_min)"
    s *= "\n  * Upper bound N_max : $(process.N_max)"
    print(io, s)
end

#--------------------------------------------------------------------------------------------------

function (process::BirthDeathProcess)(phylo::Phylogeny)
    
    lineages = copy(phylo.lineages)
    t = phylo.t
    t_end = phylo.t + process.Δt
    N = length(lineages)
    uid_count = phylo.uid_counter
    prob_division = process.λ / (process.λ + process.δ)
    
    while true
        rate = (process.λ + process.δ) * N
        t += rand(Exponential(1 / rate))
        t ≥ t_end && break
        
        n = rand(1 : N)
        if rand() < prob_division
            daughter1, daughter2, uid_count = phylo(lineages[n], uid_count, t)
            lineages[n] = daughter1
            push!(lineages, daughter2)
        else
            deleteat!(lineages, n)
        end
        
        N = length(lineages)
        process.N_min < N < process.N_max || (t_end = t; break)
    end
    
    return Phylogeny(lineages, t_end, phylo.t0, phylo.mutation_rate, uid_count)
end


function (process::BirthDeathProcess)(μ::MutationRate)
    process(μ, 1)
end

function (process::BirthDeathProcess)(μ::MutationRate, N0)
    process(μ, N0, 0.0)
end

function (process::BirthDeathProcess)(μ::MutationRate, N0, t0)
    process(Phylogeny(μ, N0; t0=t0))
end

###################################################################################################

"""
    MoranProcess{F,I} <: AbstractProcess

`MoranProcess` with turnover rate `λ` and population size `N`.

`MoranProcess` terminates after duration `Δt`.
---

## Constructors
    MoranProcess(λ::Real, N::Integer, Δt::Real)
    MoranProcess(N::Integer, Δt::Real)

Construct `MoranProcess`. If `λ` is not provided, it is set to `1.0`.
---

## Functors
    (::MoranProcess)(μ::MutationRate [, t0::Real]) -> Phylogeny
    (::MoranProcess)(μ::MutationRate, N0::Vector[, t0::Real]) -> Phylogeny
    (::MoranProcess)(phylo::Phylogeny) -> Phylogeny

Simulate `MoranProcess` and return `Phylogeny`.
"""
struct MoranProcess{F,I} <: AbstractProcess
    λ::F
    N::I
    Δt::F
    function MoranProcess(λ::Real, N::Integer, Δt::Real)
        λ, Δt, _ = promote(λ, Δt, 1.0)
        new{typeof(λ), typeof(N)}(λ, N, Δt)
    end
end

MoranProcess(N::Integer, Δt::Real) = MoranProcess(1.0, N, Δt)

function Base.show(io::IO, process::MoranProcess)
    s = "$(typeof(process))"
    s *= "\n  * Turnover rate λ   : $(process.λ)"
    s *= "\n  * Population size N : $(process.N)"
    s *= "\n  * Duration Δt       : $(process.Δt)"
    print(io, s)
end

#--------------------------------------------------------------------------------------------------

function (process::MoranProcess)(phylo::Phylogeny)
    
    lineages = copy(phylo.lineages)
    t = phylo.t
    t_stop = phylo.t + process.Δt
    uid_count = phylo.uid_counter
    rate = process.λ * process.N
    
    while t < t_stop
        t += rand(Exponential(1 / rate))
        n_divide = rand(1 : process.N)
        n_death = rand(1 : process.N)
        n_divide == n_death && continue
        
        daughter1, daughter2, uid_count = phylo(lineages[n_divide], uid_count, t)
        
        lineages[n_divide] = daughter1
        lineages[n_death] = daughter2
    end
    
    return Phylogeny(lineages, t_stop, phylo.t0, phylo.mutation_rate, uid_count)
end


function (process::MoranProcess)(μ::MutationRate, t0::Real=0.0)
    process(Phylogeny(μ, process.N; t0=t0))
end

function (process::MoranProcess)(μ::MutationRate, N0::Vector, t0::Real=0.0)
    length(N0) != process.N && error("invalid number of founder cells")
    process(Phylogeny(μ, N0; t0=t0))
end

###################################################################################################

"""
    CompositeProcess <: AbstractProcess

## Fields
* `subprocesses::Vector{AbstractProcess}` : The subprocesses.

---

## Constructor
    CompositeProcess(subprocesses...)

---

## Functor
    (::CompositeProcess)(phylo::Phylogeny) -> Phylogeny

Simulate the processes of `CompositeProcess` successively and return `Phylogeny`.
"""
struct CompositeProcess <: AbstractProcess
    subprocesses::Vector{AbstractProcess}
end

CompositeProcess(subprocesses...) = CompositeProcess([subprocesses...])

function Base.show(io::IO, process::CompositeProcess)
    s = "$(typeof(process)) with $(length(process.subprocesses)) processes"
    for (k, subprocess) in enumerate(process.subprocesses)
        s *= "\n\n$(k). $(subprocess)"
    end
    print(io, s)
end

Base.push!(process::CompositeProcess, subprocesses...) = begin
    push!(process.processes, subprocesses...)
    process
end

#--------------------------------------------------------------------------------------------------

function (process::CompositeProcess)(phylo::Phylogeny)
    for subprocess in process.subprocesses
        phylo = subprocess(phylo)
        population_size(phylo) == 0 && break
    end
    return phylo
end

###################################################################################################

# struct SelectedClone
#     process::AbstractProcess
#     t_clone::Real
#     s_clone::Real
# end

###################################################################################################