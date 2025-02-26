# QR factorization
#
# A = QR
#
# where Q is orthogonal and R is upper triangular

struct QRFactorizationSolution <: SolutionMixin
    q::AbstractMatrix
    r::AbstractMatrix
    function QRFactorizationSolution(q, r)
        sol = new(q, triu(r))

        n = size(q, 1)
        i = Matrix{Float64}(I, n, n)
        return sol
    end
end

struct QRFactorizationProblem <: ProblemMixin
    a::AbstractMatrix
end

QRFactorization = QRFactorizationProblem

function kernel(prob::QRFactorizationProblem)
    # QR factorization
    # A = QR
    #
    # where H_i is a householder reflection matrix
    # Q = H₁ H₂ ... Hₖ

    n = size(prob.a, 1)
    q = Matrix{Float64}(I, n, n)
    r = copy(prob.a)
    
    for k in 1:n    
        # Get the column vector below and including the diagonal
        v = copy(r[k:n, k])

        h = Matrix{Float64}(I, n, n)
        beta = -sign(v[1]) * norm(v)
        v[1] -= beta
        h[k:n, k:n] -= 2.0 * v * v' / dot(v, v)

        # Apply reflection to R and update Q
        r = h * r
        q = q * h
    end
    
    return QRFactorizationSolution(q, r)
end

export QRFactorization
