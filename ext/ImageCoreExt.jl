module ImageCoreExt
    using ImageCore
    using TransferFunctions
    using TransferFunctions: OpticalTransferFunction, SpatialArray
    TransferFunctions.conv(tf::OpticalTransferFunction, img::SpatialArray{T,2}) where {T<:Colorant} = mapslices(chan -> conv(tf, chan), channelview(img), dims=(2,3)) |> colorview(T)
end
