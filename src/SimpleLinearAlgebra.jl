module SimpleLinearAlgebra

    using LinearAlgebra, Printf
    using OMEinsum

    abstract type ProblemMixin end
    abstract type SolutionMixin end

    function kernel(prob::ProblemMixin)
        MethodError(kernel, (typeof(prob),)) |> throw
    end

    struct LinearAlgebraError <: Exception
        message::String
    end

    include("forward.jl")
    include("back.jl")
    include("lu.jl")
    include("qr.jl")

    export kernel
end
