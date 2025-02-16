struct ForwardSubstitutionProblem <: ProblemMixin
    l::AbstractMatrix
    b::AbstractVector
    tol::Real

    function ForwardSubstitutionProblem(l, b, tol)
        prob = new(l, b, tol)
        @assert istril(l) "matrix must be lower triangular"
        return prob
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

    for i in 1:n
        @assert abs(l[i, i]) > tol "matrix is singular"
        x[i] += b[i] / l[i, i]

        if i == 1
            continue
        end
        
        x[i] -= l[i, 1:i-1]' * x[1:i-1] / l[i, i]
    end
    return ForwardSubstitutionSolution(x)
end

export ForwardSubstitution, kernel
