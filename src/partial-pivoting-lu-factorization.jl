struct PartialPivotingLUFactorizationSolution <: SolutionMixin
    l::AbstractMatrix
    u::AbstractMatrix
    p::Vector{Int}
    function PartialPivotingLUFactorizationSolution(l, u, p)
        sol = new(l, u, p)
        @assert istril(l) "matrix must be lower triangular"
        @assert istriu(u) "matrix must be upper triangular"
        @assert isperm(p) "matrix must be a permutation matrix"
        return sol
    end
end

struct PartialPivotingLUFactorizationProblem <: ProblemMixin
    a::AbstractMatrix
    tol::Real
end

PartialPivotingLUFactorization = PartialPivotingLUFactorizationProblem

function kernel(prob::PartialPivotingLUFactorizationProblem)
    tol = prob.tol
    u = copy(prob.a)
    n = size(u, 1)

    l = zeros(n, n)
    p = collect(1:n)

    for k in 1:n-1
        v, x = findmax(i->abs(u[i, k]), k:n)
        x += k - 1

        if x != k
            u[k, :], u[x, :] = u[x, :], u[k, :]
            l[k, :], l[x, :] = l[x, :], l[k, :]
            p[k], p[x] = p[x], p[k]
        end

        if abs(u[k, k]) < tol
            continue
        end

        l[k, k] = 1.0
        l[k+1:n, k] = u[k+1:n, k] / u[k, k]
        u[k+1:n, k+1:n] -= l[k+1:n, k] * u[k, k+1:n]'
    end

    l[n, n] = 1.0
    return PartialPivotingLUFactorizationSolution(tril(l), triu(u), p)
end

export PartialPivotingLUFactorization
