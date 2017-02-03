###############################################################
######################SURVIVAL MODEL###########################
###############################################################

type SurvivalModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
    loglik::Float64
    score::Array{Float64,1}
    fischer_info::Array{Float64,2}
end


function show(io::IO, obj::SurvivalModel)
    print(io,"\nModel: ", obj.model, obj.formula,"\n\n")
    print(io,obj.coefmat)
end

###############################################################
######################EVENT####################################
###############################################################

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
######################EVENT WINDOW#############################
###############################################################

immutable EventWindow
    t₀::Float64
    t₁::Float64
end

EventWindow{S<:Real}(t::S) = EventWindow(t, t)
EventWindow(ev::Event) = ev.censored ? EventWindow(ev.time,Inf) : EventWindow(ev.time,ev.time)

function Base.show(io::IO, obj::EventWindow)
    if obj.t₀ == obj.t₁
        print(io, obj.t₀)
    elseif obj.t₁ == Inf
        print(io, obj.t₀, "+")
    else
        print(io, "[$(obj.t₀),$(obj.t₁))")
    end
end

function Base.parse(T::Type{EventWindow},str::String)
    if str[end] == '+'
        t₀ = parse(str[1:(end-1)])
        return EventWindow(t₀,Inf)
    elseif str[end] == ')'
        substr = str[2:end-1]
        splitstr = split(substr,',')
        return EventWindow(parse.(splitstr)...)
    else
        t₀ = parse(str)
        return EventWindow(t₀, t₀)
    end
end

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
        s.gradlog[i] += ders.ds_dϕ[i]
    end
    for i = 1:N
        s.gradlog[i+M] += ders.ds_dc*x[i]
    end
    for j = 1:M
        for i = 1:M
            s.heslog[i,j] += ders.d²s_dϕ²[i,j]
        end
    end
    for j = 1:M
        for i = 1:N
            s.heslog[i+M,j] += ders.d²s_dcdϕ[j]*x[i]
            s.heslog[j, i+M] = s.heslog[i+M,j]
        end
    end
    for j = 1:N
        for i = 1:N
            s.heslog[i+M,j+M] += ders.d²s_dc²*x[i]*x[j]
        end
    end
end

###############################################################
######################INTCOEFS#################################
###############################################################

type IntCoefs{R<:Number, N}
    g::Array{Vec{N,R},1}
    hggt::Array{Vec{N,R},2}
end

function IntCoefs(pdist::Distribution, degreetype = Val{50}())
    IntCoefs([clenshaw_coefs(pdist, t -> ds_dϕ(pdist,t, i), degreetype) for i in eachindex(pdist.params)],
    [clenshaw_coefs(pdist, t -> ds_dϕ(pdist,t, i)*ds_dϕ(pdist,t, j) + d²s_dϕ²(pdist,t, i, j), degreetype)
    for i in eachindex(pdist.params), j in eachindex(pdist.params)])
end
