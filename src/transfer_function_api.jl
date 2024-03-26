# TODO: Write tests <30-11-23> 
@doc raw"""
`N`-dimensional transfer function realization. You can fix the pixels sizes (`Δxy`) and sample the 
transfer function on your sensor and optical system setup.
"""
struct TransferFunctionRealization{T<:Real,N}
    "physical illumination pattern"
    transfer_function::TransferFunction{N}
    "pixel dimensions"
    Δxy::NTuple{N,Length}
end
const TFR{T,N} = TransferFunctionRealization{T,N}

function (tf::TF{N})(T::Type{<:Number}; Δxy::Union{NTuple{N,Length},Length}) where {N}
    Δxy = Δxy isa Length ? Tuple(fill(Δxy, N)) : Δxy
    TFR{T,N}(tf, Δxy)
end
(tf::TF{N})(; Δxy) where {N} = (tf)(Float64; Δxy)

# otf.jl
otf(tfr::TFR, args...; varargs...) = otf(tfr.transfer_function, args..., tfr.Δxy; varargs...)
mtf(tfr::TFR, args...; varargs...) = mtf(tfr.transfer_function, args..., tfr.Δxy; varargs...)
ptf(tfr::TFR, args...; varargs...) = ptf(tfr.transfer_function, args..., tfr.Δxy; varargs...)
otf_support(tfr::TFR, args...; varargs...) = otf_support(tfr.transfer_function, args..., tfr.Δxy; varargs...)

# psf.jl
psf(tfr::TFR, args...; varargs...) = psf(tfr.transfer_function, args..., tfr.Δxy; varargs...)
apsf(tfr::TFR, args...; varargs...) = apsf(tfr.transfer_function, args..., tfr.Δxy; varargs...)
ipsf(tfr::TFR, args...; varargs...) = ipsf(tfr.transfer_function, args..., tfr.Δxy; varargs...)

# pupil.jl ...

# TODO: Better Base.show <30-11-23> 
