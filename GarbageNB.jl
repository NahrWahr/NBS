using StaticArrays
using DiffEqBase

#####################################################################
#####################################################################
#####################################################################

abstract type Body
end

abstract type NBodySystem
end

abstract type BasicPotentialSystem <: NBodySystem
end

abstract type PotentialParameters
end

DiffEqBase.@def PosVelMass begin
    r::SVector{3, Real}
    v::SVector{3, Real}
    m::Real
end

struct MassBody <:Body
    @PosVelMass
end

struct PotentialNBodySystem{bType <: Body}
    bodies::Vector{bType}
    potentials::Dict{Symbol,<:PotentialParameters}
end

#####################################################################
#####################################################################
#####################################################################

struct GravitationalSystem{bType <: MassBody,
                           gType <: Real}<: BasicPotentialSystem
    bodies::Vector{bType}
    G::gType
end

struct GravitationalParameters{gType <: Real} <:PotentialParameters
    G::gType
end

GravitationalParameters() = GravitationalParameters(6.674e-11)

function GravAcc!(dv,
                  rs,
                  i,
                  n,
                  bodies::Vector{<:MassBody},
                  p::GravitationalParameters)

    acc = @SVector [0.0, 0.0, 0.0]
    Ri = @SVector [rs[1, i], rs[2, i], rs[3, i]]

    @inbounds for j = 1:n
        if j != 1
            Rj = @SVector [rs[1, j], rs[2, j], rs[3, j]]
            Rij = Ri - Rj
            acc -= p.G * bodies[j].m * Rij / norm(Rij)^3
        end
    end

    @. dv += acc
end

#####################################################################
#####################################################################
#####################################################################

PotentialNBodySystem(system::PotentialNBodySystem) = system

function PotentialNBodySystem(bodies::Vector{<:Body},
                              potentials::Vector{Symbol}=[])

    parameters = Dict{Symbol, PotentialParameters}()

    if :gravitational ∈ potentials
        parameters[:gravitational] = GravitationalParameters()
    end

    PotentialNBodySystem(bodies, parameters)
end

function PotentialNBodySystem(system::GravitationalSystem)

    pp = GravitationalParameters(system.G)
    potential = Dict{Symbol, PotentialParameters}(:gravitational => pp)
    PotentialNBodySystem(system.bodies, potential)
end

#####################################################################
#####################################################################
#####################################################################

abstract type BoundaryConditions
end

struct InfiniteBox{cType <: Real} <: BoundaryConditions
    boundary::SVector{6, <:cType}
end

InfiniteBox() = InfiniteBox(SVector(-Inf, Inf, -Inf, Inf, -Inf, Inf))

function GetDist(Ri,
                 Rj,
                 ::BoundaryConditions)

    Rij = Ri - Rj
    D = +(@. Rij^2...)
    return Rij, D, √D
end

#####################################################################
#####################################################################
#####################################################################

struct NBodySimulation{sType <: NBodySystem,
                       bcType <: BoundaryConditions,
                       tType <: Real,}

    system::sType
    tspan::Tuple{tType,tType}
    BndCond::bcType
    ExtGravFld
end

function NBodySimulation(system::BasicPotentialSystem,
                         tspan::Tuple{tType,tType},
                         BndCond::BoundaryConditions,
                         ExtGravFld) where {tType <: Real}

    PotSys = PotentialNBodySystem(system)
    NBodySimulation(PotSys, tspan, BndCond, ExtGravFld)
end

function NBodySimulation(system::NBodySystem,
                         tspan::Tuple{tType, tType}) where {tType <: Real}

    NBodySimulation(system, tspan, InfiniteBox(), x -> 0)
end

#####################################################################
#####################################################################
#####################################################################

function Base.show(stream::IO,
                   s::PotentialNBodySystem)
    println(stream)

    OList = [:gravitational]

    for potential in OList
        if potential ∈ keys(s.potentials)
            show(stream, s.potentials[potential])
        end
    end
end

function Base.show(stream::IO,
                   pp::GravitationalParameters)
    println(stream, "Gravitational:")
    print(stream, "\tG:"); show(stream, pp.G);
    println(stream)
end

function Base.show(stream::IO, s::NBodySimulation)
    print(stream, "Timespan: ")
    show(stream, s.tspan)
    println(stream)
    print(stream, "Boundary conditions: ")
    show(stream, s.boundary_conditions)
    println(stream)
    show(stream, s.system)
end
