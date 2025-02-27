abstract type QRFactorizationProblem <: ProblemMixin end

module QR
    using LinearAlgebra
    import SimpleLinearAlgebra.QRFactorizationProblem
    import SimpleLinearAlgebra.SolutionMixin

    struct QRFactorizationSolution <: SolutionMixin
        q::AbstractMatrix
        r::AbstractMatrix
        function QRFactorizationSolution(q, r)
            @assert istriu(r) "matrix must be upper triangular"
            @assert isapprox(q' * q, I) "matrix must be orthogonal"
            soln = new(q, r)
            return soln
        end
    end
    
    struct HouseholderReflection <: QRFactorizationProblem
        a::AbstractMatrix
        tol::Real
        function HouseholderReflection(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::HouseholderReflection)
        n = size(prob.a, 1)
        q = Matrix{Float64}(I, n, n)
        r = copy(prob.a)
        
        for i in 1:n    
            # Get the column vector below and including the diagonal
            vi = copy(r[i:n, i])

            hi = Matrix{Float64}(I, n, n)
            beta = -sign(vi[1]) * norm(vi)
            vi[1] -= beta
            hi[i:n, i:n] -= 2.0 * vi * vi' / dot(vi, vi)

            # Apply reflection to R and update Q
            r = hi * r
            q = q * hi
        end
        
        soln = QRFactorizationSolution(q, triu(r))
        return soln
    end

    struct GivenRotation <: QRFactorizationProblem
        a::AbstractMatrix
        tol::Real
        function GivenRotation(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::GivenRotation)
        n = size(prob.a, 1)
        q = Matrix{Float64}(I, n, n)
        r = copy(prob.a)

        ij = [(i, j) for i in 1:n for j in i+1:n]
        for (i, j) in ij
            x, y = r[i, i], r[j, i]
            t = atan(y, x)

            if abs(y) < prob.tol
                continue
            end

            # Manually calculate hij' * r and q * hij
            ri =  cos(t) * r[i, :] + sin(t) * r[j, :]
            rj = -sin(t) * r[i, :] + cos(t) * r[j, :]
            r[i, :], r[j, :] = ri, rj

            qi =  cos(t) * q[:, i] + sin(t) * q[:, j]
            qj = -sin(t) * q[:, i] + cos(t) * q[:, j]
            q[:, i], q[:, j] = qi, qj
        end
        return QRFactorizationSolution(q, triu(r))
    end

    struct GramSchmidt <: QRFactorizationProblem
        a::AbstractMatrix
        tol::Real
        function GramSchmidt(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::GramSchmidt)
        a = prob.a
        tol = prob.tol

        n = size(a, 1)
        q = copy(prob.a)
        r = zeros(n, n)

        r[1, 1] = norm(a[:, 1])
        @assert r[1, 1] > tol "Singular matrix"

        q[:, 1] /= r[1, 1]

        for i in 2:n
            # Method 1:
            qi = q[:, 1:i-1] # shape (n, i-1)
            ai = a[:, i] # shape (n, )
            ri = qi' * ai # shape (i-1, )

            r[1:i-1, i] = ri
            q[:, i] -= qi * ri

            r[i, i] = norm(q[:, i])
            @assert r[i, i] > tol "Singular matrix"

            q[:, i] /= r[i, i]
        end

        return QRFactorizationSolution(q, triu(r))
    end

    struct ModifiedGramSchmidt <: QRFactorizationProblem
        a::AbstractMatrix
        tol::Real
        function ModifiedGramSchmidt(a, tol=1e-10)
            prob = new(a, tol)
            return prob
        end
    end

    function kernel(prob::ModifiedGramSchmidt)
        a = copy(prob.a)
        tol = prob.tol

        n = size(a, 1)
        q = zeros(n, n)
        r = zeros(n, n)
    
        for i in 1:n
            # vi is normalized vector of a[:, i]
            # is orthogonal to all the vectors
            # before i
            ni = norm(a[:, i])
            vi = a[:, i] / ni # shape (n, )
            @assert size(vi) == (n, )

            # ai contains all the vectors
            # after i
            ai = a[:, i+1:n] # shape (n, n-i)
            # ri is the projection of vi into the 
            # basis spanned by ai
            ri = vi' * ai # shape (n-i, )
            @assert size(ri) == (1, n-i)

            r[i, i] = ni
            r[i, i+1:n] = ri
            q[:, i] = vi
            a[:, i+1:n] -= vi * ri
        end
        return QRFactorizationSolution(q, triu(r))
    end
end

kernel(prob::QRFactorizationProblem) = QR.kernel(prob)
export QR
