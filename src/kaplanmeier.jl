function kaplan_meier(events, fs, ls)
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    return getfield.(events[fs],[:time]), 1.-surv
end

function kaplan_meier(events)
    if issorted(events)
        return kaplan_meier(events, find(firsts(events)), find(lasts(events)))
    else
        sorted_events = sort(events)
        return kaplan_meier(sorted_events, find(firsts(sorted_events)), find(lasts(sorted_events)))
    end
end
