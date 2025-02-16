struct PartialPivotingLUFactorizationSolution <: SolutionMixin
    l::AbstractMatrix
    u::AbstractMatrix
    p::Vector{Int}
    function PartialPivotingLUFactorizationSolution(l, u, p)
        s = new(l, u, p)
        assert(s, istril(l), "matrix must be lower triangular")
        assert(s, istriu(u), "matrix must be upper triangular")
        assert(s, isperm(p), "matrix must be a permutation matrix")
        return s
    end
end

struct PartialPivotingLUFactorizationProblem <: ProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(p::PartialPivotingLUFactorizationProblem)
    tol = p.tol
    u = copy(p.a)
    n = size(u, 1)

    l = zeros(n, n)
    p = collect(1:n)

    for k in 1:n-1
        v, x = findmax(abs.(u[k:n, k]))
        x += k - 1

        if x != k
            # swap row k and row x of matrix u
            for i in 1:n
                u[k, i], u[x, i] = u[x, i], u[k, i]
            end

            # swap row k and row x of matrix l
            for i in 1:n
                l[k, i], l[x, i] = l[x, i], l[k, i]
            end

            # swap row k and row x of matrix p
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

export PartialPivotingLUFactorizationProblem
