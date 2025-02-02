module SimpleLinearAlgebra

    using LinearAlgebra

    abstract type LinearSystemProblemMixin end
    abstract type LinearSystemSolutionMixin end

    function kernel(prob::LinearSystemProblemMixin)
        MethodError(kernel, (typeof(prob),)) |> throw
    end

    include("forward_substitution.jl")
    include("back_substitution.jl")
    include("lufact.jl")
end
