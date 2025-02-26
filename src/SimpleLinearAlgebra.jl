module SimpleLinearAlgebra

    using LinearAlgebra, Printf
    import Printf.@sprintf  # Fixed macro import

    abstract type ProblemMixin end
    abstract type SolutionMixin end

    function kernel(prob::ProblemMixin)
        MethodError(kernel, (typeof(prob),)) |> throw
    end

    struct LinearAlgebraError <: Exception
        message::String
    end

    include("forward-substitution.jl")
    include("back-substitution.jl")
    include("lu-factorization.jl")
    include("partial-pivoting-lu-factorization.jl")
    include("qr-factorization.jl")
end
