type EventHistoryModel
    model::AbstractString
    formula::Formula
    #eventtype::DataType
    #coef::Vector{Float64}
    coefmat::CoefTable
    M::ModelFrame
    #IC::Matrix{Float64}
    #invhess::Matrix{Float64}
    #vcov::Matrix{Float64}
    #opt::Vector
    #grad::Vector{Float64}
    #X::Matrix{Float64}
    #eventtime::Matrix
    #chaz::Matrix{Float64}
end

# coefnames !!!

function show(io::IO, obj::EventHistoryModel)
    print(io,"\nModel: ", obj.model, "; ", obj.formula,"\n")
    #n = size(obj.eventtime,1)
    #events::Int = sum(obj.eventtime[:,2])
    #print(io,"\nn=",n,", events=",events,"\n\n")
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


type xyfunc{S<:Real, T<:Real}
    x::Vector{S}
    y::Vector{T}
end

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










#Altra funzione importante somme parziali!
# Avanti e indietro per somme parziali!!!!

# function forward_sum!{T}(v, times, censored, addend::AbstractArray{T})
#     for i = 1:length(times)
#         if censored[i]
#             v[i,:] = zero(T)
#         elseif (i==1) || (times[i] > times[i-1])
#             v[i,:] = addend[i,:]
#         else
#             v[i,:] = v[i-1,:]+addend[i,:]
#         end
#     end
#     return
# end
#
# function backward_sum!{T}(w, v::AbstractArray{T}, times, censored)
#     for i = length(times):-1:1
#         if censored[i]
#             w[i,:] = zero(T)
#         elseif (i==size(times,1)) || (times[i+1] > times[i]) || censored[i+1]
#             w[i,:] = v[i,:]
#         else
#             w[i,:] = w[i+1,:]
#         end
#     end
#     return
# end
#
# function total_partial_sum!(w,v,times,censored,addend)
#     forward_sum!(v,times,censored,addend)
#     backward_sum!(w,v,times,censored)
#     return
# end
#
# #Ordina DataFrame!! # da ottimizzare ma serve questa funzioneeee!
# function extra_info(times, censored)
#     dd = DataFrame()
#     dd[:time] = times
#     dd[:censored] = censored
#     # Order by event_time, putting censored after normals
#     sort!(dd; cols = [:time, :censored])
#     # Add number of currently alive processes
#     dd[:alive] = collect(size(dd,1):-1:1)
#
#     #Count Death Number!
#     dd[:death_num] = zeros(Int64, size(dd,1))
#     # Check number of simulataneous deaths
#     dd[:death_tot] = zeros(Int64, size(dd,1))
#     total_partial_sum!(dd[:death_tot],dd[:death_num], times, censored, ones(Int64, size(dd,1)) )
#     return dd
# end
#
#
# function add_prob!(km)
#     km[:surv] = zeros(size(km,1))
#     surv = 1.
#     for i in 1:size(km,1)
#         surv *= (1-km[i,:death_tot]/km[i,:alive])
#         km[i,:surv] = surv
#     end
#     km[:prob] = -diff(vcat(1.,km[:surv]))
#     return
# end
#
#
# function kaplan_meier(times::AbstractArray{S,1}, censored::AbstractArray{Bool,1})
#     dd = extra_info(times, censored)
#     km = dd[dd[:death_num].== 1, :]
#     add_prob!(km)
#     return km
# end
