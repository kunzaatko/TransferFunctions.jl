# TODO: Add support for defining interfaces on types that are not from the `TransferFunctions` base module. See the
# `Estimation` module <05-05-25> 
# TODO: Add support for adding documentation to the interface <25-04-25> 
# TODO: Support giving an output type to the interface function i.e. f(a::Abstract, b::ConcreteIn)::ConcreteOut should work <24-04-25> 
# TODO: Support only indicating the interface type i.e. f(::Abstract, b::ConcreteIn) should work <24-04-25> 
# TODO: Support parametric types i.e. f(::Abstract{A}, b::ConcreteIn{B}) should work <24-04-25> 
macro require_interface(ex)
    if ex.head == :(=) || ex.head == Symbol("function")
        signature = ex.args[1]
    else
        signature = ex
    end

    (typeof(signature) == Expr && signature.head == :call) || throw(ArgumentError("`@require_interface` may be used only on a function or a `:call` expression."))
    interface_name = signature.args[1]
    length(signature.args) > 1 || throw(ArgumentError("In `@require_interface`, an interface must have atleast one argument."))
    interface_params = signature.args[2:end]
    (typeof(first(interface_params)) == Expr && first(interface_params).head == :(::)) || throw(ArgumentError("In `@require_interface`, the first argument of the interface must have a known type otherwise the interfacing type is unknown."))
    interfacing_type = interface_params[1].args[2]
    isabstracttype(eval(interfacing_type)) || throw(ArgumentError("In `@require_interface`, the interfacing type must be abstract."))
    rest_params = length(signature.args) > 1 ? interface_params[2:end] : nothing

    return esc(quote
        function $interface_name(A::T, $(rest_params...)) where {T<:$interfacing_type}
            error("`$(nameof(T))` does not implement the obligatory interface `$($interface_name)($(join($interface_params, ", ")))`!")
        end
    end)
end

include("types/extension-interface.jl")
include("types/transfer-function-interface.jl")
