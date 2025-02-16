struct LUFactorizationSolution <: SolutionMixin
    l::AbstractMatrix
    u::AbstractMatrix

    function LUFactorizationSolution(l, u)
        s = new(l, u)
        @assert istril(l) "matrix must be lower triangular"
        @assert istriu(u) "matrix must be upper triangular"
        return s
    end
end

abstract type LUFactorizationProblemMixin <: ProblemMixin end

struct Version1 <: LUFactorizationProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(prob::Version1)
    tol = prob.tol
    u = copy(prob.a)
    n = size(u, 1)

    l = Matrix{Float64}(I, n, n)

    for k in 1:n-1
        @assert abs(u[k, k]) > tol "Gaussian elimination failed"

        m_k = Matrix{Float64}(I, n, n)
        for i in k+1:n
            m_k[i, k] = -u[i, k] / u[k, k]
        end

        l = l * inv(m_k)
        u = m_k * u
    end

    return LUFactorizationSolution(tril(l), triu(u))
end

struct Version2 <: LUFactorizationProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(prob::Version2)
    tol = prob.tol
    u = copy(prob.a)
    n = size(u, 1)

    # initialize l to be the identity matrix
    l = Matrix{Float64}(I, n, n)

    for k in 1:n-1
        @assert abs(u[k, k]) > tol "Gaussian elimination failed"
        l[k+1:n, k] = u[k+1:n, k] / u[k, k]
        u[k+1:n, k+1:n] -= l[k+1:n, k] * u[k, k+1:n]'
    end

    return LUFactorizationSolution(tril(l), triu(u))
end

LUFactorizationV1 = Version1
LUFactorizationV2 = Version2
export LUFactorizationV1, LUFactorizationV2

LUFactorization = LUFactorizationV2
export LUFactorization
