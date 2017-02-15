###############################################################
######################SURVIVAL MODEL###########################
###############################################################

abstract SurvivalModel

function show(io::IO, obj::SurvivalModel)
    print(io,"\nModel: ", obj.model, obj.formula,"\n\n")
    print(io,obj.coefmat)
end

coef(SM::SurvivalModel) = SM.coefmat.cols[1]

immutable CoxModel <: SurvivalModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
    loglik::Float64
    score::Array{Float64,1}
    fischer_info::Array{Float64,2}
end

immutable AftModel{D<:Distribution} <: SurvivalModel
    model::AbstractString
    formula::Formula
    coefmat::CoefTable
    M::ModelFrame
    loglik::Float64
    score::Array{Float64,1}
    fischer_info::Array{Float64,2}
    dist::D
end

coef_reg(AM::AftModel) = coef(AM)[length(AM.dist.params)+1:end]
coef_dist(AM::AftModel) = AM.dist.params

function predict(AM::AftModel)
    model_matrix = DataFrames.ModelMatrix(AM.M)
    X = model_matrix.m
    coefs = coef_reg(AM)
    exp.(X*coefs)
end

function rand(AM::AftModel)
    preds = predict(AM)
    [rand(AM.dist)*pred for pred in preds]
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
######################AFTRESP##################################
###############################################################
immutable AftResp{R<:Number, D<:Distribution}
    W::Array{R,2}
    X::Array{R,2}
    Xβ::Array{R,1}
    Θ::Array{R,1}
    dist::D
    M::Int64
    N::Int64
    I1::Array{Int64,1}
    I2::Array{Int64,1}
    τs::Array{R,2}
    cdfs::Array{R,2}
    Δcdf::Array{R,1}
    loglik::Array{R,1}
end

function AftResp(S, X, dist)
    M = length(dist.params)
    N = size(X,2)
    Xβ = zeros(size(X,1))
    Θ = ones(size(X,1))
    L = length(S)
    W = [(k == 1 ? s.t₀ : s.t₁) for s in S, k in 1:2]
    instants = [s.t₁ == s.t₀ for s in S]
    I1 = find(instants)
    I2 = find(!instants)
   return AftResp(W, X, Xβ, Θ, dist, M, N, I1, I2,
                   zeros(L,2), ones(L,2),
                   zeros(L), zeros(L))
end
###############################################################
######################INTCOEFS#################################
###############################################################

immutable IntCoefs{R<:Number, N}
    aux::Array{R,1}
    g::Array{Vec{N,R},1}
    hggt::Array{Vec{N,R},2}
end

function IntCoefs{N}(pdist::Distribution, degreetype::Val{N} = Val{50}())
    aux = auxvec(pdist)
    int_coefs = IntCoefs(aux,[Vec{N, Float64}(zeros(N)) for i in eachindex(pdist.params)],
    [Vec{N, Float64}(zeros(N)) for i in eachindex(pdist.params), j in eachindex(pdist.params)])
    for i in eachindex(pdist.params)
        int_coefs.g[i] = clenshaw_coefs(pdist, t -> ds(pdist,t, int_coefs, i), degreetype)
    end
    for i in eachindex(pdist.params), j in eachindex(pdist.params)
        int_coefs.hggt[i, j] = clenshaw_coefs(pdist, t -> ds(pdist,t, int_coefs, i)*ds(pdist,t,int_coefs, j) +
        d²s(pdist,t, int_coefs,i, j), degreetype)
    end
    return int_coefs
end
