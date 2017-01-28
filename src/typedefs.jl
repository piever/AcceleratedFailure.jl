###############################################################
######################EVENT###################################
###############################################################

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

###############################################################
######################SURVWINDOW###############################
###############################################################

type SurvWindow{S<:Real, T<:Real}
    t1::S
    t2::T
    instant::Bool
end

SurvWindow{S<:Real}(t::S) = SurvWindow(t,true)
SurvWindow{S<:Real}(t::S,i::Bool) = i ? SurvWindow(t,t,true) : SurvWindow(t,Inf,false)
SurvWindow{S<:Real, T<:Real}(s::S,t::T) = SurvWindow(s,t,false)

function Base.show(io::IO, obj::SurvWindow)
    print(io, obj.instant ? obj.t1:"[$(obj.t1),$(obj.t2))")
end

type Pdistribution{F1<:Function, F2<:Function, F3<:Function, F4<:Function, F5<:Function}
    pdf::F1
    cdf::F2
    M::Int64
    gradlog::F3
    heslog::F4
    quantile::F5
end

Pdistribution(a,b,c) = Pdistribution(a,b,c,
(ϕ,t)->error("undefined gradlog!"), (ϕ,t)->error("undefined gradhes!"), (ϕ,t)->error("undefined quantile!"))

ll(s::SurvWindow,pdist::Pdistribution,ϕ) = s.instant ?
log(pdist.pdf(ϕ,s.t1)) : log(pdist.cdf(ϕ,s.t2)-pdist.cdf(ϕ,s.t1))

function ll(s::SurvWindow,pdist::Pdistribution,ϕ, coef)
    s.instant ? ll(coef*s,pdist,ϕ) + log(coef) : ll(coef*s,pdist,ϕ)
end

Base.:*(a,b::SurvWindow) = SurvWindow(a*b.t1, a*b.t2, b.instant)
