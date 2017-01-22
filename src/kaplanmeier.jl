function kaplan_meier(events, fs, ls)
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    return xyfunc(getfield.(events[fs],[:time]), 1.-surv)
end

function kaplan_meier(events, sorted::Bool = false)
    sorted ? kaplan_meier(events, find(firsts(events)), find(lasts(events))) :
             kaplan_meier(sort(events),true)
end
