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

        for i in 1:n
            for j in i+1:n
                x, y = r[i, i], r[j, i]
                c = x / sqrt(x^2 + y^2)
                s = y / sqrt(x^2 + y^2)

                if abs(s) < prob.tol
                    continue
                end

                h = Matrix{Float64}(I, n, n)
                h[[i, j], [i, j]] = [c -s; s c]

                r = h' * r
                q = q * h
            end
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
            # q[:, i] = a[:, i]
            
            for j in 1:i-1
                r[j, i] = dot(q[:, j], a[:, i])
                q[:, i] -= q[:, j] * r[j, i]
            end

            r[i, i] = norm(q[:, i])
            @assert r[i, i] > tol "Singular matrix"

            q[:, i] /= r[i, i]
        end

        return QRFactorizationSolution(q, triu(r))
    end
end

kernel(prob::QRFactorizationProblem) = QR.kernel(prob)
export QR
