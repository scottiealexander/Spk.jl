module SpkTuning

using SpkCore, Statistics, Optim#,NLopt,

include("./TuningUtils.jl")

export TuningData, fit, f0, f1, interpolate
export Contrast, LogisticContrast, Area, SSE, WeightedSSE

# ============================================================================ #
f0(x::AbstractVector{<:Real}, tf::Real=0.0) = mean(x)
f1(x::AbstractVector{<:Real}, tf::Real) = abs(f1xfm(cycle_mean(x, floor(Int, 1.0 / tf))))
# ============================================================================ #
function trial_reduce(op::Function, data::TuningData, bin_size::Real=0.001, calc_error::Bool=false)
    kbin = 1:round(Int, data.dur / bin_size)
    ulab = sort(unique(data.lab))
    mn = zeros(length(ulab))
    err = ones(length(ulab))
    for k in eachindex(ulab)
        kevt = findall(isequal(ulab[k]), data.lab)
        if calc_error
            d = psth(data.ts, data.evt[kevt], kbin, bin_size)[1] ./ bin_size
            tmp = zeros(size(d, 2))
            kc = 1
            for col in eachcol(tmp)
                tmp[kc] = op(col, data.tf * bin_size)
                kc += 1
            end
            # tmp = mapslices(x->op(x, data.tf), d, dims=2)
            mn[k] = mean(tmp)
            err[k] = std(tmp)
        else
            d = vec(mean(psth(data.ts, data.evt[kevt], kbin, bin_size)[1], dims=2)) ./ bin_size
            mn[k] = op(d, data.tf * bin_size)
        end
    end
    return FitData(ulab, mn, err)
end
# ============================================================================ #
function fit(::Type{T}, data::TuningData, op::Function, bin_size::Real=0.001) where {T<:TuningFunction{SSE}}
    return fit(T, trial_reduce(op, data, bin_size))
end
# ---------------------------------------------------------------------------- #
function fit(::Type{T}, data::TuningData, op::Function, bin_size::Real=0.001) where {T<:TuningFunction{WeightedSSE}}
    return fit(T, trial_reduce(op, data, bin_size, true))
end
# ---------------------------------------------------------------------------- #
function fit(::Type{T}, fd::FitData) where {T<:TuningFunction}
    p0, lb, ub = parameters(T, fd)
    count = 0
    fobj(p::Vector) = begin
        count += 1
        return evaluate(T, fd, p)
    end

    # Optim.jl version using LBFGS
    res = optimize(fobj, lb, ub, p0, Fminbox(LBFGS()); autodiff = :forward)
    flag = Optim.x_converged(res) || Optim.f_converged(res) || Optim.g_converged(res) ?
        :converged : :converge_fail
    return TuningResult(T, fd, Optim.minimizer(res), Optim.minimum(res), flag)

    # opt = Opt(:LN_PRAXIS, length(p0))
    #
    # lower_bounds!(opt, lb)
    # upper_bounds!(opt, ub)
    # min_objective!(opt, fobj)
    # minf, minx, ret = optimize(opt, p0)
    #
    # return TuningResult(T, fd, minx, minf, ret)
end
# ============================================================================ #
function objective(::Type{T}, fd::FitData, model::Vector) where {T<:TuningFunction{SSE}}
    return sum((model - fd.y).^2)
end
# ---------------------------------------------------------------------------- #
function objective(::Type{T}, fd::FitData, model::Vector) where {T<:TuningFunction{WeightedSSE}}
    return sum(((model - fd.y) .^ 2.0) ./ fd.err)
end
# ---------------------------------------------------------------------------- #
function evaluate(::Type{T}, fd::FitData, p::Vector) where {T<:TuningFunction}
    return objective(T, fd, generate(T, fd.x, p))
end
# ---------------------------------------------------------------------------- #
function interpolate(tr::TuningResult{T}, n::Int=100, log::Bool=false) where {T<:TuningFunction}

    xmn, xmx = extrema(tr.data.x)

    if log
        x = 10 .^ range(log10(xmn), log10(xmx), length=n)
    else
        x = range(xmn, xmx, length=n)
    end
    y = generate(T, collect(x), tr.p)

    return x, y
end
# ============================================================================ #
# function test()
#     db = get_database("contrast", id->id>100);
#     ret, lgn, evt, lab = get_data(db, 1);
#     dur = get_uniform_param(db[1], "stimulus_duration");
#     tf = get_uniform_param(db[1], "temporal_frequency");
#     td = TuningData(ret, evt, lab, dur, tf);
#     tr = fit(Contrast{SSE}, td, f0, 0.001);
#     x, y = interpolate(tr, 100);
#     plot(tr.data.x, tr.data.y, ".")
#     plot(x, y)
# end
# ============================================================================ #
end
