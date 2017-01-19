# what's wrong with the types?

function nelson_aalen(res::EventHistoryModel)
    nelson_aalen(collect(res.M.df[:,1]),find(firsts(res.M.df[:,1])),
    find(lasts(res.M.df[:,1])),(DataFrames.ModelMatrix(res.M).m)[:,2:end],res.coefmat.cols[1])
end

function nelson_aalen(events, fs, ls, X, β)
    Xβ = X*β
    Θ = exp.(Xβ)
    ds = ls-fs+1
    ns = after(Θ)[fs]
    chaz = cumsum(ds./ns)
    return xyfunc(getfield.(events[fs],[:time]), chaz)
end

function nelson_aalen(events, fs, ls)
    ds = ls-fs+1
    ns = length(events)- fs +1
    chaz = cumsum(ds./ns)
    return xyfunc(getfield.(events[fs],[:time]), chaz)
end


function nelson_aalen(events, sorted::Bool = false)
    sorted ? nelson_aalen(events, find(firsts(events)), find(lasts(events))) :
             nelson_aalen(sort(events),true)
end

function nelson_aalen{S<:Real, T<:Real}(dist1::AbstractVector{S},dist2::AbstractVector{T})
    times = min.(dist1,dist2)
    censored = (dist1 .> dist2)
    return nelson_aalen(Event.(times,censored))
end

function nelson_aalen(tdist1 :: Distributions.Distribution, tdist2 :: Distributions.Distribution, n :: Int64)
    dist1,dist2 = rand.([tdist1,tdist2],[n])
    return nelson_aalen(dist1,dist2)
end
