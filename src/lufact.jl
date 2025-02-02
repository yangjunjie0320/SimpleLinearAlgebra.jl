struct LUFactorizationSolution <: SolutionMixin
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

struct LUFactVersion1 <: ProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(p::LUFactVersion1)
    tol = p.tol
    u = copy(p.a)
    n = size(u, 1)

    l = zero(u)
    l[1:n+1:end] .= 1

    for k in 1:n-1
        if abs(u[k, k]) < tol
            ArgumentError("Gaussian elimination failed") |> throw
        end

        m_k = zero(u)
        m_k[1:n+1:end] .= 1
        for i in k+1:n
            m_k[i, k] = -u[i, k] / u[k, k]
        end

        l = l * inv(m_k)
        u .= m_k * u
    end

    return LUFactorizationSolution(tril(l), triu(u))
end



struct LUFactVersion2 <: ProblemMixin
    a::AbstractMatrix
    tol::Real
end

function kernel(p::LUFactVersion2)
    tol = p.tol
    u = copy(p.a)
    n = size(u, 1)

    l = zero(u)
    l[1:n+1:end] .+= 1

    for k in 1:n-1
        if abs(u[k, k]) < tol
            ArgumentError("Gaussian elimination failed") |> throw
        end

        for i in k+1:n
            l[i, k] = u[i, k] / u[k, k]
        end

        for j in k+1:n
            for i=k+1:n
                u[i, j] -= l[i, k] * u[k, j]
            end
        end
    end

    return LUFactorizationSolution(tril(l), triu(u))
end

export LUFactVersion1, LUFactVersion2

LUFactorizationProblem = LUFactVersion2
export LUFactorizationProblem
