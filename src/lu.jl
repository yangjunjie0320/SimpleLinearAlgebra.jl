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

    struct GaussianEliminationV1 <: LUFactorizationProblemMixin
        a::AbstractMatrix
        tol::Real
        function GaussianEliminationV1(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::GaussianEliminationV1)
        tol = prob.tol
        u = copy(prob.a)
        n = size(u, 1)

        l = Matrix{Float64}(I, n, n)
        # l * u = a now, we insert resolution of the identity matrix
        # l * inv(lk) * lk * u to make l * inv(lk) lower triangular,
        # and lk * u to make it upper triangular
        # as lk is constructed from gaussian elimination, lk * lu
        # will eliminate the k-th column of u
        # and as inv(lk) is also lower triangular, the product of
        # l * inv(lk) is lower triangular

        for i in 1:n-1
            uii = u[i, i]
            @assert abs(uii) > tol "Gaussian elimination failed"
            
            # lk is lower triangular matrix
            li = Matrix{Float64}(I, n, n)
            li[i+1:n, i] -= u[i+1:n, i] / uii
            li_inv = inv(li)

            l = l * li_inv
            u = li * u
        end

        soln = LUFactorizationSolution(tril(l), triu(u))
        return soln
    end

    struct GaussianEliminationV2 <: LUFactorizationProblemMixin
        a::AbstractMatrix
        tol::Real
        function GaussianEliminationV2(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::GaussianEliminationV2)
        tol = prob.tol
        u = copy(prob.a)
        n = size(u, 1)

        # initialize l to be the identity matrix
        l = Matrix{Float64}(I, n, n)

        for i in 1:n-1
            uii = u[i, i]
            @assert abs(uii) > tol "Gaussian elimination failed"

            ci = u[i+1:n, i] / uii # shape (n-i, ), the i-th column of U
            ri = u[i, i+1:n] # shape (1, n-i), the i-th row of U
            
            l[i+1:n, i] = ci
            u[i+1:n, i+1:n] -= ci * ri'
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
        
        # initialize l to be the identity matrix
        l = zeros(n, n)
        p = collect(1:n)
    
        for i in 1:n-1
            # v is the max element, x is the index
            v, j = findmax(u[i:n, i])
            j += i - 1
    
            if i != j # if the max element is not on the diagonal
                u[i, :], u[j, :] = u[j, :], u[i, :]
                l[i, :], l[j, :] = l[j, :], l[i, :]
                p[i], p[j] = p[j], p[i]
            end
            
            uii = u[i, i]
            if abs(uii) < tol
                continue
            end
            
            l[i, i] = 1.0
            ri = u[i, i+1:n]
            ci = u[i+1:n, i] / uii
            l[i+1:n, i] = ci
            u[i+1:n, i+1:n] -= ci * ri'
        end
        
        l[n, n] = 1.0
        soln = PartialPivotingLUFactorizationSolution(tril(l), triu(u), p)
        return soln
    end
end

kernel(prob::LUFactorizationProblemMixin) = LU.kernel(prob)

LUFactorization = LU.GaussianEliminationV2
export LUFactorization

struct CholeskyFactorizationProblem <: ProblemMixin
    a::AbstractMatrix
    tol::Real
    function CholeskyFactorizationProblem(a, tol=1e-10)
        prob = new(a, tol)
        return prob
    end
end

struct CholeskyFactorizationSolution <: SolutionMixin
    l::AbstractMatrix
    function CholeskyFactorizationSolution(l)
        @assert istril(l) "matrix must be lower triangular"
        soln = new(l)
        return soln
    end
end

function kernel(prob::CholeskyFactorizationProblem)
    tol = prob.tol
    l = copy(prob.a)
    n = size(l, 1)

    for i in 1:n
        @assert l[i, i] > tol "Cholesky factorization failed"
        l[i, i] = sqrt(l[i, i])
        l[i+1:n, i] /= l[i, i]

        li = l[i+1:n, i]
        l[i+1:n, i+1:n] -= li * li'
    end

    soln = CholeskyFactorizationSolution(tril(l))
    return soln
end

CholeskyFactorization = CholeskyFactorizationProblem
CholeskyDecomposition = CholeskyFactorization
export CholeskyFactorization, CholeskyDecomposition
