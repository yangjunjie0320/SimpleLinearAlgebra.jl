struct LUFactorizationProblem <: LinearSystemProblemMixin
    a::AbstractMatrix
    tol::Real
end

struct LUFactorizationSolution <: LinearSystemSolutionMixin
    l::AbstractMatrix
    u::AbstractMatrix

    function LUFactorizationSolution(l, u)
        if !istril(l)
            ArgumentError("Matrix must be lower triangular") |> throw
        end

        if !istriu(u)
            ArgumentError("Matrix must be upper triangular") |> throw
        end

        return new(l, u)
    end
end

LUFactorization = LUFactorizationProblem

function elemination_step(a::AbstractMatrix{T}, k::Integer) where T
    # given a matrix A, perform the k-th elimination step.
    # M A = N, which M is a lower triangular matrix with ones on the diagonal,
    # and N is the matrix A with the k-th column eliminated.
    # The inverse of M is simply given by Minv = 2I - m

    n = size(a, 1)
    @assert size(a, 2) == n
    @assert 1 <= k <= n-1

    m = identity_matrix(T, n)
    for i in k+1:n
        m[i, k] = -a[i, k] / a[k, k]
    end
    return m
end

function kernel(prob::LUFactorizationProblem)
    u = copy(prob.a)
    n = size(u, 1)
    @assert size(u, 2) == n

    tol = prob.tol

    l = zero(u)

    l[1:n+1:end] .= 1

    for k in 1:n-1

        if abs(u[k, k]) < tol
            ArgumentError("Gaussian elimination failed") |> throw
        end

        for i=k+1:n
            l[i, k] = u[i, k] / u[k, k]
        end

        for j in k+1:n
            for i in k+1:n
                u[i, j] -= l[i, k] * u[k, j]
            end
        end
    end

    return LUFactorizationSolution(l, triu(u))
end

export LUFactorization, kernel