module SimpleLinearAlgebra

    using LinearAlgebra

    abstract type LinearSystemProblemMixin end
    abstract type LinearSystemSolutionMixin end

    function kernel(prob::LinearSystemProblemMixin)
        throw(MethodError(kernel, (prob,)))
    end

    include("forward_substitution.jl")
    include("back_substitution.jl")
end
