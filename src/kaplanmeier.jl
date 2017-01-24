function kaplan_meier(events::Array, fs, ls)
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    return getfield.(events[fs],[:time]), 1.-surv
end

function kaplan_meier(events::Array)
    if issorted(events)
        return kaplan_meier(events, find(firsts(events)), find(lasts(events)))
    else
        sorted_events = sort(events)
        return kaplan_meier(sorted_events, find(firsts(sorted_events)), find(lasts(sorted_events)))
    end
end

kaplan_meier(events::NullableArray, args...) = kaplan_meier(events.values[!events.isnull], args...)
kaplan_meier(events::AbstractArray, args...) = kaplan_meier(convert(Array,events), args...)
