struct LUFactorizationSolution <: SolutionMixin
    l::AbstractMatrix
    u::AbstractMatrix

    function LUFactorizationSolution(l, u)
        s = new(l, u)
        assert(s, istril(l), "matrix must be lower triangular")
        assert(s, istriu(u), "matrix must be upper triangular")
        return s
    end
end

abstract type LUFactorizationProblemMixin <: ProblemMixin end

struct Version1 <: LUFactorizationProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(p::Version1)
    tol = p.tol
    u = copy(p.a)
    n = size(u, 1)

    l = Matrix{Float64}(I, n, n)

    for k in 1:n-1
        assert(p, abs(u[k, k]) > tol, "Gaussian elimination failed")

        m_k = Matrix{Float64}(I, n, n)
        for i in k+1:n
            m_k[i, k] = -u[i, k] / u[k, k]
        end

        l = l * inv(m_k)
        u .= m_k * u
    end

    return LUFactorizationSolution(tril(l), triu(u))
end

struct Version2 <: LUFactorizationProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(p::Version2)
    tol = p.tol
    u = copy(p.a)
    n = size(u, 1)

    # initialize l to be the identity matrix
    l = Matrix{Float64}(I, n, n)

    for k in 1:n-1
        assert(p, abs(u[k, k]) > tol, "Gaussian elimination failed")
        l[k+1:n, k] = u[k+1:n, k] / u[k, k]
        u[k+1:n, k+1:n] -= l[k+1:n, k] * u[k, k+1:n]'
    end

    return LUFactorizationSolution(tril(l), triu(u))
end

LUFactVersion1 = Version1
LUFactVersion2 = Version2
export LUFactVersion1, LUFactVersion2

LUFactorizationProblem = Version2
export LUFactorizationProblem
