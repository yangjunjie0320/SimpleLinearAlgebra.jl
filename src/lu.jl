abstract type LUFactorizationProblemMixin <: ProblemMixin end

module LU
    using LinearAlgebra
    import SimpleLinearAlgebra.LUFactorizationProblemMixin
    import SimpleLinearAlgebra.SolutionMixin

    struct LUFactorizationSolution <: SolutionMixin
        l::AbstractMatrix
        u::AbstractMatrix
        function LUFactorizationSolution(l, u, p=nothing)
            @assert istril(l) "matrix must be lower triangular"
            @assert istriu(u) "matrix must be upper triangular"
            soln = new(l, u)
            return soln
        end
    end

    struct Version1 <: LUFactorizationProblemMixin
        a::AbstractMatrix
        tol::Real
        function Version1(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::Version1)
        tol = prob.tol
        u = copy(prob.a)
        n = size(u, 1)

        l = Matrix{Float64}(I, n, n)

        for k in 1:n-1
            @assert abs(u[k, k]) > tol "Gaussian elimination failed"

            m_k = Matrix{Float64}(I, n, n)
            m_k[k+1:n, k] = -u[k+1:n, k] / u[k, k]

            l = l * inv(m_k)
            u = m_k * u
        end

        soln = LUFactorizationSolution(tril(l), triu(u))
        return soln
    end

    struct Version2 <: LUFactorizationProblemMixin
        a::AbstractMatrix
        tol::Real
        function Version2(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
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

        soln = LUFactorizationSolution(tril(l), triu(u))
        return soln
    end

    struct PartialPivotingLUFactorizationSolution <: SolutionMixin
        l::AbstractMatrix
        u::AbstractMatrix
        p::AbstractVector
        function PartialPivotingLUFactorizationSolution(l, u, p)
            @assert istril(l) "matrix must be lower triangular"
            @assert istriu(u) "matrix must be upper triangular"
            @assert isperm(p) "permutation vector must be a permutation"
            soln = new(l, u, p)
            return soln
        end
    end

    struct PartialPivoting <: LUFactorizationProblemMixin
        a::AbstractMatrix
        tol::Real
        function PartialPivoting(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::PartialPivoting)
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
        soln = PartialPivotingLUFactorizationSolution(tril(l), triu(u), p)
        return soln
    end
end

kernel(prob::LUFactorizationProblemMixin) = LU.kernel(prob)

LUFactorization = LU.Version2
export LUFactorization
