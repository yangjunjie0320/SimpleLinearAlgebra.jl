struct BackSubstitutionProblem <: LinearSystemProblemMixin
    u::AbstractMatrix 
    b::AbstractVector
    tol::Real

    function BackSubstitutionProblem(u, b, tol)
        if !istriu(u)
            ArgumentError("Matrix must be upper triangular") |> throw
        end
        return new(u, b, tol)
    end
end

struct BackSubstitutionSolution <: LinearSystemSolutionMixin
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
    x[n] = b[n] / u[n, n]

    for i in n-1:-1:1
        if abs(u[i, i]) < tol
            throw(ArgumentError("Matrix is singular"))
        end
        x[i] = (b[i] - sum(u[i, j] * x[j] for j in i+1:n)) / u[i, i]
    end

    return BackSubstitutionSolution(x)
end


export BackSubstitution, kernel