################################################################################
## Type 'EventClass' for Event History objects
## Subtypes:
##   Surv:      (Right)-censored survival outcome
##   SurvInt:   Interval-censored outcome (Left,Right,Interval)
##   SurvTrunc: Right-censoring + left-truncation
##   CompRisk:  Competing risks model
##
## constructor: Event
## Access methods: Time, Status, Cause, ...
################################################################################

###{{{ Types/classes: Surv,CompRisk,...

abstract EventClass
#abstract EventInt <: EventClass

immutable Surv <: EventClass
    Time::Number  # Exit time
    Status::Bool  # Censoring status
end



# immutable SurvInt <: EventClass
#     Time::Number     # Time 1
#     Time2::Number    # Interval [Time;Time2]. Use Inf/-Inf for left/right censoring
#     Status::Int # Censoring (0:none,1:right,2:left,3:interval)
#     function SurvInt(Time,Time2)
#         if Time>Time2 error("time 1 larger than time 2") end
#         if (Time==Time2)
#             status = EventHistory.CensNot
#         elseif Time2==Inf
#             status = EventHistory.CensRight
#         elseif Time==-Inf
#             status = EventHistory.CensLeft
#         else
#             status = EventHistory.CensInt
#         end
#         new(Time,Time2,status)
#     end
# end

###}}} Types

###{{{ show methods
function show(io::IO, obj::Surv)
    print(io, obj.Time, obj.Status>0 ? "":"+")
end
function show(io::IO, obj::SurvTrunc)
    print(io, "(", obj.Entry, ";", obj.Time, obj.Status>0 ? "":"+","]")
end
