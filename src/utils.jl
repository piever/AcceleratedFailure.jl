type Event{T<:Real}
    time::T
    censored::Bool
end

function Base.isless(a::Event, b::Event)
    if a.time < b.time
        return true
    elseif a.time > b.time
        return false
    else
        return isless(a.censored, b.censored)
    end
end

#Altra funzione importante somme parziali!
# Avanti e indietro per somme parziali!!!!

function forward_sum!{T}(v, times, censored, addend::AbstractArray{T})
    for i = 1:length(times)
        if censored[i]
            v[i,:] = zero(T)
        elseif (i==1) || (times[i] > times[i-1])
            v[i,:] = addend[i,:]
        else
            v[i,:] = v[i-1,:]+addend[i,:]
        end
    end
    return
end

function backward_sum!{T}(w, v::AbstractArray{T}, times, censored)
    for i = length(times):-1:1
        if censored[i]
            w[i,:] = zero(T)
        elseif (i==size(times,1)) || (times[i+1] > times[i]) || censored[i+1]
            w[i,:] = v[i,:]
        else
            w[i,:] = w[i+1,:]
        end
    end
    return
end

function total_partial_sum!(w,v,times,censored,addend)
    forward_sum!(v,times,censored,addend)
    backward_sum!(w,v,times,censored)
    return
end

#Ordina DataFrame!! # da ottimizzare ma serve questa funzioneeee!
function extra_info(times, censored)
    dd = DataFrame()
    dd[:time] = times
    dd[:censored] = censored
    # Order by event_time, putting censored after normals
    sort!(dd; cols = [:time, :censored])
    # Add number of currently alive processes
    dd[:alive] = collect(size(dd,1):-1:1)

    #Count Death Number!
    dd[:death_num] = zeros(Int64, size(dd,1))
    # Check number of simulataneous deaths
    dd[:death_tot] = zeros(Int64, size(dd,1))
    total_partial_sum!(dd[:death_tot],dd[:death_num], times, censored, ones(Int64, size(dd,1)) )
    return dd
end
