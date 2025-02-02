struct ForwardSubstitutionProblem <: ProblemMixin
    l::AbstractMatrix
    b::AbstractVector
    tol::Real

    function ForwardSubstitutionProblem(l, b, tol)
        if !istril(l)
            ArgumentError("Matrix must be lower triangular") |> throw
        end
        return new(l, b, tol)
    end
end

ForwardSubstitution = ForwardSubstitutionProblem

struct ForwardSubstitutionSolution <: SolutionMixin
    x::AbstractVector
end

function kernel(prob::ForwardSubstitutionProblem)
    l = prob.l
    b = copy(prob.b)
    n = length(b)
    tol = prob.tol

    x = zeros(n)
    x[1] = b[1] / l[1, 1]
    for i in 2:n
        if abs(l[i, i]) < tol
            ArgumentError("Matrix is singular") |> throw
        end
        x[i] = (b[i] - sum(l[i, j] * x[j] for j in 1:i-1)) / l[i, i]
    end
    return ForwardSubstitutionSolution(x)
end

export ForwardSubstitution, kernel
