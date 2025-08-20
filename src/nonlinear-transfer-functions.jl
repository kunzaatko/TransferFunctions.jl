"""
    NonLinearTransferFunction{N} <: TransferFunction{N}
A supertype for all non-linear transfer functions.

See also [`TransferFunction`](@ref), [`LinearShiftInvariantTransferFunction`](@ref)
"""
abstract type NonLinearTransferFunction{N} <: TransferFunction{N} end
