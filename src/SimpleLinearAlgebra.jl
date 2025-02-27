module SimpleLinearAlgebra
    using LinearAlgebra, Printf

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
    include("eig.jl")

    export kernel
end
