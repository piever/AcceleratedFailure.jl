function kaplan_meier(sort_ev::SortedEvents)
    fs, ls, events = sort_ev.fs, sort_ev.ls, sort_ev.events
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    v = surv.^2.*cumsum(ds./(ns.*(ns-ds)))
    return getfield.(events[fs],[:time]), surv, v
end

kaplan_meier(events::Array, sorted = false) = kaplan_meier(SortedEvents(events, sorted))

if isdefined(:NullableArray)
    kaplan_meier(events::NullableArray, args...) = kaplan_meier(events.values[!events.isnull], args...)
end
kaplan_meier(events::AbstractArray, args...) = kaplan_meier(convert(Array,events), args...)
