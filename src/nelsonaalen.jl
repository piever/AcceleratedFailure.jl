function nelson_aalen(events, fs, ls, Θ = ones(length(events)))
    ds = ls-fs+1
    ns = after(Θ)[fs]
    chaz = cumsum(ds./ns)
    return getfield.(events[fs],[:time]), chaz
end

function nelson_aalen(events, Θ = ones(length(events)))
    if issorted(events)
        return nelson_aalen(events, find(firsts(events)), find(lasts(events)),Θ)
    else
        p = sortperm(events)
        sorted_events = events[p]
        sorted_Θ = Θ[p]
        return nelson_aalen(sorted_events, find(firsts(sorted_events)), find(lasts(sorted_events)), sorted_Θ)
    end
end

function nelson_aalen(res::EventHistoryModel)
    nelson_aalen(convert(Array,res.M.df[:,1]),
    exp.(((DataFrames.ModelMatrix(res.M).m)[:,2:end])*(res.coefmat.cols[1])))
end
