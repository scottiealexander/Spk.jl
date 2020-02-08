
export filter_events_by_duration

# ============================================================================ #
function filter_events_by_duration(evt::Vector{<:Real},
            dur::AbstractFloat, nevt::Integer, pad::AbstractFloat=0.05)

    if length(evt) < nevt
        error("Number of event timestamps is < the given number of event")
    end

    pad = 1.0 .+ [-pad; +pad]
    df = diff(evt)
    bgood = (df .>= dur*pad[1]) .& (df .<= dur*pad[2])
    ngood = sum(bgood)

    if ngood == nevt-1
        evt_ts = evt[[bgood; true]]

    elseif ngood == nevt
        evt_ts = evt[[bgood; false]]

    else
        push!(bgood, true)

        # proportion of events of duration dur
        pdur = ngood ./ length(df)

        if pdur <= .1
            # few / no events are of length dur: events *MUST* only mark stim
            # onset
            total_dur = nanmedian(df)
            buse = (df .>= total_dur*pad[1]) & (df .<= total_dur*pad[2])
            nuse = sum(buse)

            if nuse == nevt
                buse = [buse; false]
            elseif nuse == nevt-1
                buse = [buse; true]
            else
                error("Failed to locate the correct number of events")
            end

            evt_ts = evt[buse]

        elseif pdur >= .9
            if ngood + 1 >= nevt*2
                # most / all events are of length dur: events mark on and offset
                # of stim and blank_dur == stim_dur
                evt_ts = evt[bgood]

                #index of first stim on event (by definition precedes stim off)
                buse = falses(length(evt_ts))
                buse[1:2:(nevt*2)] .= true
                evt_ts = evt_ts[buse]
            elseif ngood > nevt
                # events mark stim on but we have extra events of approx. the
                # correct duration so just take the first nevt...
                @warn("Correct event could not be unambiguously identified...")
                @warn("falling back to heuristic method. Check the results!")
                evt_ts = evt[bgood]
                evt_ts = evt_ts[1:nevt]
            else
                error("Failed to locate the correct number of events")
            end
        else
            error("Event durations are too inconsistent to estimate blank duration")
        end
    end
    return evt_ts
end
# ============================================================================ #
nanmedian(x::Vector{<:AbstractFloat}) = median(x[!isnan(x)])
# ============================================================================= #
