type EventHistoryModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
    loglik::Float64
    score::Array{Float64,1}
    fischer_info::Array{Float64,2}
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

function Base.parse(T::Type{Event},str::String)
    if str[end] == '+'
        return Event(parse(str[1:(end-1)]), true)
    else
        return Event(parse(str), false)
    end
end

Base.isless(a::Event, b::Event) = isless((a.time, a.censored), (b.time,b.censored))

# Compute first and last! Could be optimized!
firsts(S) = [!S[t].censored && (t==1 || S[t] > S[t-1]) for t = 1:length(S)]
lasts(S) = [!S[t].censored && (t==length(S) || S[t+1] > S[t]) for t = 1:length(S)]


# Computes "complementary cumulative sum" on first dimension
function after(v)
    afterv = init_after(v)
    after!(afterv,v)
    return afterv
end


function init_after(v)
    newsize = collect(size(v))
    newsize[1] += 1
    afterv = zeros(eltype(v),newsize...)
    return afterv
end

function after!{N,T}(afterv, v::AbstractArray{T,N})
    trails = fill(:,N-1)
    cumsum!(@view(afterv[(end-1):-1:1,trails...]),@view(v[end:-1:1,trails...]),1)
    return
end
