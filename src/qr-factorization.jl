struct QRFactorizationSolution <: SolutionMixin
    q::AbstractMatrix
    r::AbstractMatrix
    function QRFactorizationSolution(q, r)
        sol = new(q, r)
        @assert isortho(q) "matrix must be orthogonal"
        @assert istriu(r) "matrix must be upper triangular"
        return sol
    end
end

struct QRFactorizationProblem <: ProblemMixin
    a::AbstractMatrix
end

QRFactorization = QRFactorizationProblem

function householder_reflection!(q, r)
    m, n = size(r)
    
    if m == 1
        return q, r
    else
        h = Matrix{Float64}(I, m, m)
        for i in 1:m-1
            h[i+1:m, i] = -2 * r[i+1:m, i] / r[i, i]
        end
        q = h * q
        r = h * r
        householder_reflection!(q, r)
    end
end

function kernel(prob::QRFactorizationProblem)
    r = copy(prob.a)
    n = size(r, 1)
    q = Matrix{Float64}(I, n, n)    


end

export PartialPivotingLUFactorization
