abstract type LinearTransferFunction <: TransferFunction end
conv(t::LinearTransferFunction, ::SpatialImage) = no_implemementation_error(typeof(t), :conv)
deconv(t::LinearTransferFunction, ::SpatialImage) = no_implemementation_error(typeof(t), :deconv)
transfer(t::LinearTransferFunction, img::SpatialImage) = conv(t, img)
restore(t::LinearTransferFunction, img::SpatialImage) = deconv(t, img)

include("otf/optical-transfer-function.jl")
include("psf/point-spread-function.jl")
