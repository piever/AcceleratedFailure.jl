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
    print(io,"\nModel: ", obj.model, obj.formula,"\n\n")
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
    t₀::S
    t₁::T
    instant::Bool
end

SurvWindow{S<:Real}(t::S) = SurvWindow(t,true)
SurvWindow{S<:Real}(t::S,i::Bool) = i ? SurvWindow(t,t,true) : SurvWindow(t,Inf,false)
SurvWindow{S<:Real, T<:Real}(s::S,t::T) = SurvWindow(s,t,false)

function Base.show(io::IO, obj::SurvWindow)
    print(io, obj.instant ? obj.t₀:"[$(obj.t₀),$(obj.t₁))")
end

###############################################################
######################PDISTRIBUTIONS###########################
###############################################################

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
log(pdist.pdf(ϕ,s.t₀)) : log(pdist.cdf(ϕ,s.t₁)-pdist.cdf(ϕ,s.t₀))

function ll(s::SurvWindow,pdist::Pdistribution,ϕ, coef)
    s.instant ? ll(coef*s,pdist,ϕ) + log(coef) : ll(coef*s,pdist,ϕ)
end

Base.:*(a,b::SurvWindow) = SurvWindow(a*b.t₀, a*b.t₁, b.instant)

###############################################################
######################DERIVATIVES AFT##########################
###############################################################

type Derivatives{R<:Number}
    loglik::R
    ds_dϕ::Array{R,1}
    ds_dc::R
    d²s_dϕ²::Array{R,2}
    d²s_dc²::R
    d²s_dcdϕ::Array{R,1}
end

Derivatives(M::Int64) = Derivatives(0., zeros(M), 0., zeros(M,M), 0., zeros(M))

type SmoothLog{R<:Number}
    val::R
    gradlog::Array{R,1}
    heslog::Array{R,2}
end

function zero!{R<:Number}(s::SmoothLog{R})
    s.val = zero(R)
    s.gradlog[:] = zero(R)
    s.heslog[:] = zero(R)
    return
end

function add_ders!(s::SmoothLog, ders::Derivatives, x, M , N)
    s.val += ders.loglik
    for i = 1:M
        s.grad[i] += ders.ds_dϕ[i]
    end
    for i = 1:N
        s.grad[i+M] += ders.ds_dc*x[i]
    end
    for j = 1:M
        for i = 1:M
            s.hes[i,j] += ders.d²s_dϕ²[i,j]
        end
    end
    for j = 1:M
        for i = 1:N
            s.hes[i+M,j] += ders.d²s_dcdϕ[i,j]
            s.hes[j, i+M] = s.hes[i+M,j]
        end
    end
    for j = 1:N
        for i = 1:N
            s.hes[i+M,j+M] += ders.d²s_dc²[i,j]*x[i]*x[j]
        end
    end
end
