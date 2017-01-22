# what's wrong with the types?

function nelson_aalen(res::EventHistoryModel)
    nelson_aalen(collect(res.M.df[:,1]),find(firsts(res.M.df[:,1])),
    find(lasts(res.M.df[:,1])),(DataFrames.ModelMatrix(res.M).m)[:,2:end],res.coefmat.cols[1])
end

nelson_aalen(events, fs, ls, X, β) = nelson_aalen(events, fs, ls, exp.(X*β))

function nelson_aalen(events, fs, ls, Θ = ones(length(events)))
    ds = ls-fs+1
    ns = after(Θ)[fs]
    chaz = cumsum(ds./ns)
    return xyfunc(getfield.(events[fs],[:time]), chaz)
end

function nelson_aalen(events, sorted::Bool = false)
    sorted ? nelson_aalen(events, find(firsts(events)), find(lasts(events))) :
             nelson_aalen(sort(events),true)
end
