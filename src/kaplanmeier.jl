function kaplan_meier{R}(events::AbstractVector{Event{R}}, fs, ls)
    n = length(events)
    ds = ls-fs+1
    ns = length(events)- fs +1
    surv = cumprod(1.-ds./ns)
    return xyfunc(getfield.(events[fs],[:time]), 1.-surv)
end


function kaplan_meier{R}(events::AbstractVector{Event{R}}, sorted::Bool = false)
    sorted ? kaplan_meier(events, find(firsts(events)), find(lasts(events))) :
             kaplan_meier(sort(events),true)
end

function kaplan_meier{S<:Real, T<:Real}(dist1::AbstractVector{S},dist2::AbstractVector{T})
    times = min.(dist1,dist2)
    censored = (dist1 .> dist2)
    return kaplan_meier(Event.(times,censored))
end

function kaplan_meier(tdist1 :: Distributions.Distribution, tdist2 :: Distributions.Distribution, n :: Int64)
    dist1,dist2 = rand.([tdist1,tdist2],[n])
    return kaplan_meier(dist1,dist2)
end



# Compute first and last! Could be optimized!
firsts{R}(S::Vector{Event{R}}) = [!S[t].censored && (t==1 || S[t] > S[t-1]) for t = 1:length(S)]
lasts{R}(S::Vector{Event{R}}) = [!S[t].censored && (t==length(S) || S[t+1] > S[t]) for t = 1:length(S)]

#f = find(firsts)
#l = find(lasts)

















# get_mean(km; event_time = :event_time) = (km[event_time]' *km[:prob])[1]
# get_std(km; event_time = :event_time) = sqrt(((km[event_time].^2)' * km[:prob])[1]-
# get_mean(km; event_time = event_time)^2)
# function get_both(km; event_time = :event_time)
#     media = get_mean(km; event_time = :event_time)
#     standard = sqrt(((km[event_time].^2)' * km[:prob])[1]-media^2)
#     return media, standard
# end
# function get_ratio(args...;kwargs...)
#     m, σ = get_both(args...;kwargs...)
#     return σ/m
# end
