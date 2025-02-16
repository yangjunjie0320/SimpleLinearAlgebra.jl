struct BackSubstitutionProblem <: ProblemMixin
    u::AbstractMatrix 
    b::AbstractVector
    tol::Real

    function BackSubstitutionProblem(u, b, tol)
        prob = new(u, b, tol)
        assert(prob, istriu(u), "matrix must be upper triangular")
        return prob
    end
end

struct BackSubstitutionSolution <: SolutionMixin
    x::AbstractVector
end

BackSubstitution = BackSubstitutionProblem

function kernel(prob::BackSubstitutionProblem)
    u = prob.u
    b = copy(prob.b)
    n = length(b)
    tol = prob.tol

    # solve the system
    x = zeros(n)

    for i in n:-1:1
        assert(prob, abs(u[i, i]) > tol, "matrix is singular")
        x[i] += b[i] / u[i, i]

        if i == n
            continue
        end

        x[i] -= u[i, i+1:n]' * x[i+1:n] / u[i, i]
    end

    return BackSubstitutionSolution(x)
end

export BackSubstitution, kernel
