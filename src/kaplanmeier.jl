function add_prob!(km)
    km[:surv] = zeros(size(km,1))
    surv = 1.
    for i in 1:size(km,1)
        surv *= (1-km[i,:death_tot]/km[i,:alive])
        km[i,:surv] = surv
    end
    km[:prob] = -diff(vcat(1.,km[:surv]))
    return
end

function kaplan_meier{S<:Real}(times::AbstractArray{S,1}, censored::AbstractArray{Bool,1})
    dd = extra_info(times, censored)
    km = dd[dd[:death_num].== 1, :]
    add_prob!(km)
    return km
end

function kaplan_meier{S<:Real, T<:Real}(dist1::AbstractArray{S,1},dist2::AbstractArray{T,1})
    event = min.(dist1,dist2)
    censored = (dist1 .> dist2)
    return kaplan_meier(event,censored)
end

function kaplan_meier(tdist1 :: Distributions.Distribution, tdist2 :: Distributions.Distribution, n :: Int64)
    dist1,dist2 = rand.([tdist1,tdist2],[n])
    return kaplan_meier(dist1,dist2)
end



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
