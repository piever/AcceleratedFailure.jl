function nelson_aalen(sort_ev::SortedEvents, Θ = ones(length(sort_ev.events)))
    fs, ls, events = sort_ev.fs, sort_ev.ls, sort_ev.events
    ds = ls-fs+1
    ns = after(Θ)[fs]
    chaz = cumsum(ds./ns)
    v = cumsum(ds./ns.^2)
    return [event.time for event in events[fs]], chaz, v
end

nelson_aalen(events::Array, Θ = ones(length(events))) =
nelson_aalen(SortedEvents(events), Θ)

if isdefined(:NullableArray)
    nelson_aalen(events::NullableArray, args...) = nelson_aalen(events.values[!events.isnull], args...)
end
nelson_aalen(events::AbstractArray, args...) = nelson_aalen(convert(Array,events), args...)

function nelson_aalen(res::CoxModel)
    nelson_aalen(convert(Array,res.M.df[:,1]),
    exp.(((DataFrames.ModelMatrix(res.M).m)[:,2:end])*(res.coefmat.cols[1])))
end

function survival(res::CoxModel)
    x_chaz,y_chaz,v_chaz = nelson_aalen(res::CoxModel)
    x, y, v = x_chaz, exp.(-y_chaz), exp.(-2y_chaz) .* v_chaz
    return x, y, v
end

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
