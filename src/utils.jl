type EventHistoryModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
end

# coefnames !!!

function show(io::IO, obj::EventHistoryModel)
    print(io,"\nModel: ", obj.model, "; ", obj.formula,"\n\n")
    print(io,obj.coefmat)
end

immutable Event{T<:Real}
    time::T
    censored::Bool
end

Event(x) = Event(x,false)

function Base.show(io::IO, obj::Event)
    print(io, obj.time, obj.censored ? "+":"")
end

Base.isless(a::Event, b::Event) = isless((a.time, a.censored), (b.time,b.censored))

# Compute first and last! Could be optimized!
firsts(S) = [!S[t].censored && (t==1 || S[t] > S[t-1]) for t = 1:length(S)]
lasts(S) = [!S[t].censored && (t==length(S) || S[t+1] > S[t]) for t = 1:length(S)]

# Preallocate??
function after{N,T}(v::AbstractArray{T,N})
    cv = cumsum(v,1)
    newsize = collect(size(cv))
    newsize[1] += 1
    afterv = zeros(eltype(cv),newsize...)
    if N==1
        afterv[1:end-1] = cv[end] - cv + v
    elseif N==2
        afterv[1:end-1,:] = (cv[end:end,:] .- cv) + v
    elseif N==3
        afterv[1:end-1,:,:] = (cv[end:end,:,:] .- cv) + v
    else
        error("$v can be at most 3-dimensional!")
    end
    return afterv
end
