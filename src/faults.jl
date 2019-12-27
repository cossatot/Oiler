module Faults


using Parameters

using ..Oiler: VelocityVectorSphere

@with_kw struct Fault
    trace::Float64

end