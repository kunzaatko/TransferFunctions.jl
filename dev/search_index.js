var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = TransferFunctions","category":"page"},{"location":"#TransferFunctions","page":"Home","title":"TransferFunctions","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for TransferFunctions.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [TransferFunctions]","category":"page"},{"location":"#TransferFunctions.TransferFunctions","page":"Home","title":"TransferFunctions.TransferFunctions","text":"This module implements diffraction transfer function models that are widely used in microscopy. The base type for any transfer function is TransferFunction and further disambiguates to MeasuredTransferFunction and ModelTransferFunction.\n\nwarning: Warning\nAll of the models assume incoherent illumination sources (i.e. where the radiation from the source is inocoherent). This is an approximation for most optical setups, but one that is used almost universally (e.g. in fluorescence microscopy).\n\n\n\n\n\n","category":"module"},{"location":"#TransferFunctions.BornWolf","page":"Home","title":"TransferFunctions.BornWolf","text":"Born & Wolf model of the transfer function for a circular aperture.\n\nParameters\n\n* `λ::Unitful.Length`: (emission) wavelength \n* `nᵢ::Number`: index of refraction of the immersion medium\n* `NA::Number`: numerical aperture\n\nExtended help\n\nThe Born & Wolf model is a scalar diffraction model derived for perfect systems. It assumes that the only aberration of the system is due to defocus. Modern microscope objectives are designed to provide optimal imaging conditions for sources located directly on the coverslip, in which case the Born & Wolf model is applicable (if the coverslip and immersion is used as designed). The model disregards spherical and higher order aberrations that are due to the source of illumination being shifted from the coverslip boundary.\n\n\n\n\n\n","category":"type"},{"location":"#TransferFunctions.GibsonLanni","page":"Home","title":"TransferFunctions.GibsonLanni","text":"Gibson & Lanni model of the transfer function for a circular aperture.\n\nParameters\n\nExtended help\n\nThe Gibson & Lanni model is assumes that, disregarding defocus, all observed aberrations are generated by factors external to the objective (i.e. originating in the sample, coverslip and immersion medium combination). These aberrations can be characterized by the optical path difference between a ray in a perfect system (see BornWolf) and a ray under experimental conditions.\n\n\n\n\n\n","category":"type"},{"location":"#TransferFunctions.SymmetricPupilFunction","page":"Home","title":"TransferFunctions.SymmetricPupilFunction","text":"If the pupil function of the system is symmetric, the OTF as well as the PSF are radially symmteric which can be used to optimize the calculations\n\n\n\n\n\n","category":"type"},{"location":"#TransferFunctions.apsf","page":"Home","title":"TransferFunctions.apsf","text":"Amplitude point spread function\n\n\n\n\n\n","category":"function"},{"location":"#TransferFunctions.apsf-Tuple{TransferFunctions.TransferFunction, Vararg{Any}}","page":"Home","title":"TransferFunctions.apsf","text":"apsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}\napsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.ipsf","page":"Home","title":"TransferFunctions.ipsf","text":"Intensity point spread function\n\n\n\n\n\n","category":"function"},{"location":"#TransferFunctions.ipsf-Tuple{TransferFunctions.TransferFunction, Vararg{Any}}","page":"Home","title":"TransferFunctions.ipsf","text":"ipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix{<:Real}\nipsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix{<:Real}\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.mtf-Tuple","page":"Home","title":"TransferFunctions.mtf","text":"modulation transfer function\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.otf","page":"Home","title":"TransferFunctions.otf","text":"optical transfer function\n\n\n\n\n\n","category":"function"},{"location":"#TransferFunctions.otf-Tuple{TransferFunctions.ClosedFormOTFModel, Tuple{Integer, Integer}, Tuple{Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}}}","page":"Home","title":"TransferFunctions.otf","text":"otf(tf, wh, [Δxy]; δ=(0, 0))\notf(tf, img, [Δxy];...)\n\nGenerate an otf for the given transfer function with the size wh (size of img) and with a pixel distance of Δxy\n\nnote: Note\nThe pixel size/distance (Δxy) is only required for a model transfer function. If it is  not provided for a MeasuredTransferFunction, then the pixel distance from the measurement is used.\n\nArguments\n\ntf::TransferFunction: transfer function model/measure to generate the OTF for\nwh::Tuple{Integer, Integer} or wh::Integer: (width, height) of the generated OTF. wh ↦ (wh, wh) if wh isa Integer.\nΔxy::Tuple{Length, Length} or Δxy::Length: Separation of pixels in the x and y dimensions of the generated   OTF image. Δxy ↦ (Δxy, Δxy) if Δxy isa Length.\nδ::Tuple = (0,0): shift of the OTF in the image plane in pixels. This is useful for some algorithms, e.g. in    structured illumination microscopy reconstruction algorithms.\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.psf","page":"Home","title":"TransferFunctions.psf","text":"Intensity point spread function i.e. the intensity ratio and phase shift of the sample intensity density\n\n\n\n\n\n","category":"function"},{"location":"#TransferFunctions.psf-Union{Tuple{TF}, Tuple{Type{TransferFunctions.SymmetricPupilFunction{TF}}, TF, Tuple{Integer, Integer}, Tuple{Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}}}} where TF<:TransferFunctions.ClosedFormPSFModel","page":"Home","title":"TransferFunctions.psf","text":"psf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Tuple{Length,Length})::OffsetMatrix\npsf(tf::TransferFunction, wh::Tuple{Integer,Integer}, Δxy::Length)::OffsetMatrix\n\nGenerate the psf with size wh with Δxy being the distance between the samples in the x and y dimensions\n\njulia> tf = BornWolf(488u\"nm\", 1.7, 1.7)\nBornWolf{Float64}(488.0 nm, 1.7, 1.7)\n\ntf = BornWolf(488u\"nm\", 1.7, 1.7) # hide\njulia> psf(tf, (11,11), 60u\"nm\")\njulia> psf(tf, (5,5), (40u\"nm\", 50u\"nm\")) # different pixelsizes in x and y direction\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.psf-Union{Tuple{TF}, Tuple{Type{TransferFunctions.SymmetricPupilFunction{TF}}, TF, Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}, Union{Quantity{T, 𝐋, U}, Level{L, S, Quantity{T, 𝐋, U}} where {L, S}} where {T, U}}} where TF<:TransferFunctions.ClosedFormPSFModel","page":"Home","title":"TransferFunctions.psf","text":"psf(tf::TransferFunction, x::Length, y::Length)\n\nSample the PSF of the transfer function model/data at the point (x,y)\n\njulia> tf = BornWolf(488u\"nm\", 1.4, 1.7)\nBornWolf{Float64}(488.0 nm, 1.4, 1.7)\n\ntf = BornWolf(488u\"nm\", 1.4, 1.7) # hide\npsf(tf, 0u\"nm\", 5u\"nm\")\npsf.(tf, 0u\"nm\", -400u\"nm\":100u\"nm\":400u\"nm\")\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.ptf-Tuple","page":"Home","title":"TransferFunctions.ptf","text":"phase transfer function\n\n\n\n\n\n","category":"method"},{"location":"#TransferFunctions.pupil-Tuple{TransferFunctions.TransferFunction}","page":"Home","title":"TransferFunctions.pupil","text":"pupil function\n\n\n\n\n\n","category":"method"}]
}
