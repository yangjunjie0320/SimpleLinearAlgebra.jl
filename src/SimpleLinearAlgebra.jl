module SimpleLinearAlgebra

    using LinearAlgebra

    abstract type ProblemMixin end
    abstract type SolutionMixin end

    function kernel(prob::ProblemMixin)
        MethodError(kernel, (typeof(prob),)) |> throw
    end

    include("forward_substitution.jl")
    include("back_substitution.jl")
    include("lufact.jl")
end
