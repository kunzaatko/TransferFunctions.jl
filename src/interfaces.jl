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

    @assert signature.head == :call "`@require_interface` may be used only on a function or a `:call` expression"
    interface_name = signature.args[1]
    @assert length(signature.args) > 1 "An interface must have atleast one argument"
    interface_params = signature.args[2:end]
    @assert first(interface_params).head == :(::) "First argument of the interface must have a known type otherwise the interfacing type is unknown"
    interfacing_type = interface_params[1].args[2]
    @assert isabstracttype(eval(interfacing_type)) "The interfacing type must be abstract"
    rest_params = length(signature.args) > 1 ? interface_params[2:end] : nothing

    return esc(quote
        function $interface_name(A::T, $(rest_params...)) where {T<:$interfacing_type}
            error("`$(nameof(T))` does not implement the obligatory interface `$($interface_name)($(join($interface_params, ", ")))`!")
        end
    end)
end

include("types/extension-interface.jl")
include("types/transfer-function-interface.jl")
