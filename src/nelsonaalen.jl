function nelson_aalen(events::Array, fs, ls, Θ = ones(length(events)))
    ds = ls-fs+1
    ns = after(Θ)[fs]
    chaz = cumsum(ds./ns)
    return getfield.(events[fs],[:time]), chaz
end

function nelson_aalen(events::Array, Θ = ones(length(events)))
    if issorted(events)
        return nelson_aalen(events, find(firsts(events)), find(lasts(events)),Θ)
    else
        p = sortperm(events)
        sorted_events = events[p]
        sorted_Θ = Θ[p]
        return nelson_aalen(sorted_events, find(firsts(sorted_events)), find(lasts(sorted_events)), sorted_Θ)
    end
end

nelson_aalen(events::NullableArray, args...) = nelson_aalen(events.values[!events.isnull], args...)
nelson_aalen(events::AbstractArray, args...) = nelson_aalen(convert(Array,events), args...)

function nelson_aalen(res::EventHistoryModel)
    nelson_aalen(convert(Array,res.M.df[:,1]),
    exp.(((DataFrames.ModelMatrix(res.M).m)[:,2:end])*(res.coefmat.cols[1])))
end

chaz2cdf(x,y) = (x,1. -exp.(-y))

chaz2haz(x,y) = (x, diff(vcat(0.,y)))

function chaz2haz(x,y,bw; npoints = 1000)
    pts = linspace(0., maximum(x), npoints)
    ycont = zeros(size(pts))
    index = vcat(0,round.([Int64],x/step(pts)))+1
    yext = vcat(0.,y)
    for i in 1:length(x)
        ycont[index[i]:index[i+1]] = yext[i]
    end
    dist = Normal(0.,bw)
    smoother = pdf(dist, quantile(dist,0.001):step(pts):quantile(dist,0.999))
    append!(ycont,fill(y[end],div(length(smoother),2)+1))
    ysmooth = conv(ycont,smoother)[div(length(smoother),2):(length(pts)+div(length(smoother),2))]
    smoothaz = diff(ysmooth)
    return pts, smoothaz
end
