function kaplan_meier(sort_ev::SortedEvents)
    fs, ls, events = sort_ev.fs, sort_ev.ls, sort_ev.events
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    return getfield.(events[fs],[:time]), surv
end

kaplan_meier(events::Array) = kaplan_meier(SortedEvents(events))

if isdefined(:NullableArray)
    kaplan_meier(events::NullableArray, args...) = kaplan_meier(events.values[!events.isnull], args...)
end
kaplan_meier(events::AbstractArray, args...) = kaplan_meier(convert(Array,events), args...)
