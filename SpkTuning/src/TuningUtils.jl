# ============================================================================ #
abstract type AbstractObjective end
abstract type SSE <: AbstractObjective end
abstract type WeightedSSE <: AbstractObjective end

abstract type TuningFunction{T<:AbstractObjective} end
# ============================================================================ #
mutable struct TuningData
    ts::Vector{Float64}
    evt::Vector{Float64}
    lab::Vector{Float64}
    dur::Float64
    tf::Float64
end
# ============================================================================ #
mutable struct FitData
    x::Vector{Float64}
    y::Vector{Float64}
    err::Vector{Float64}
end
# ---------------------------------------------------------------------------- #
FitData(x::Vector{Float64}, y::Vector{Float64}) = FitData(x, y, zeros(length(x)))
# ============================================================================ #
mutable struct TuningResult{T}
    res::Dict{Symbol, Float64}
    p::Vector{Float64}
    min::Float64
    flag::Symbol
    data::FitData
end
# ---------------------------------------------------------------------------- #
function TuningResult(::Type{T}, fd::FitData, p::Vector{Float64},
    mn::Float64, flag::Symbol) where {T<:TuningFunction}

    return TuningResult{T}(Dict(k=>v for (k,v) in zip(labels(T), p)), p, mn, flag, fd)
end
# ---------------------------------------------------------------------------- #
Base.getindex(tr::TuningResult, s::Symbol) = tr.res[s]
# ============================================================================ #
struct Contrast{T} <: TuningFunction{T} end
labels(::Type{Contrast{T}}) where T = [:max, :n, :c50, :baseline]
function parameters(::Type{Contrast{T}}, fd::FitData) where T
    xmn, xmx = extrema(fd.x)
    ymn, ymx = extrema(fd.y)

    p0 = [ymx, 1.5, xmn + ((xmx - xmn)*0.3), ymn]

    lb = [ymn, 0.6, xmn, 0.0]
    ub = [ymx*1.2, 6.0, xmx, ymn * 1.2]

    return p0, lb, ub
end
function generate(::Type{Contrast{T}}, x::Vector{Float64}, p::Vector) where T
    return p[1] .* (x.^p[2] ./ (x.^p[2] .+ p[3].^p[2])) .+ p[4]
end
# ============================================================================ #
struct LogisticContrast{T} <: TuningFunction{T} end
labels(::Type{LogisticContrast{T}}) where T = [:max, :k, :c50, :baseline]
function parameters(::Type{LogisticContrast{T}}, fd::FitData) where T
    p0, lb, ub = parameters(Contrast{SSE}, fd)
    p0[2] = 1.0
    # ub[2] = 15.0
    return p0, lb, ub
end
function generate(::Type{LogisticContrast{T}}, x::Vector{Float64}, p::Vector) where T
    return ((p[1] - p[4]) ./ (1.0 .+ exp.(-p[2].*(x .- p[3])))) .+ p[4]
end
# ============================================================================ #
struct Area{T} <: TuningFunction{T} end
labels(::Type{Area{T}}) where T = [:ke, :ki, :se, :si]
function parameters(::Type{Area{T}}, fd::FitData) where T
    mx = maximum(fd.y)
    p0 = [maximum(fd.y), 1.0, 1.0, 4.0]
    lb = [0.0, 0.0, fd.x[1], fd.x[1]]
    ub = [Inf, Inf, fd.x[end], fd.x[end]]
    return p0, lb, ub
end
function generate(::Type{Area{T}}, x::Vector{Float64}, p::Vector) where T
    step = 0.01
    ge = zeros(length(x))
    gi = copy(ge)
    for k = 1:length(x)
        x_sub = 2.0 .* range(-x[k]/2, x[k]/2, step=step)
        ge[k] = trapz(exp.(-1 .*((x_sub) ./ p[3]).^2))
        gi[k] = trapz(exp.(-1 .*((x_sub) ./ p[4]).^2))
    end
    return p[1].*ge.*step .- p[2].*gi.*step
end
# ============================================================================ #
trapz(y::AbstractVector{<:Real}) = 0.5 * y[1] + sum(y[2:end-1]) + 0.5 * y[end]
# ============================================================================ #
