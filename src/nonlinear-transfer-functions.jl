"""
    NonLinearShiftInvariantTransferFunction <: TransferFunction
A supertype for all non-linear transfer functions.

See also [`TransferFunctions`](@ref), [`LinearShiftInvariantTransferFunction`](@ref)
"""
abstract type NonLinearTransferFunction <: TransferFunction end

# TODO: Create a structure for this type branch <01-08-25> 
"""
    ImpulseResponseMapping <: NonLinearTransferFunction
An impulse response mapping is a description of a transfer function that prescribes an impulse response to every point
in the object plane. This may be a [`point spread function`](@ref PointSpreadFunction) or some other function that
prescribes the response to the object plane point in the image plane.

An transfer defined by an impulse response mapping is a linear system, but is not translation invariant hence does not
belong to the class of linear transfer functions.
"""
struct ImpulseResponseMapping <: NonLinearTransferFunction
end
