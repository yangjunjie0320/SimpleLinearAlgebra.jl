struct BackSubstitutionProblem <: ProblemMixin
    u::AbstractMatrix 
    b::AbstractVector
    tol::Real

    function BackSubstitutionProblem(u, b, tol=1e-10)
        prob = new(u, b, tol)
        @assert istriu(u) "matrix must be upper triangular"
        return prob
    end
end

struct BackSubstitutionSolution <: SolutionMixin
    x::AbstractVector
end

function kernel(prob::BackSubstitutionProblem)
    u = prob.u
    b = copy(prob.b)
    n = length(b)
    tol = prob.tol

    # solve the system
    x = zeros(n)

    for i in n:-1:1
        uii = u[i, i]
        @assert abs(uii) > tol "matrix is singular"
        x[i] += b[i] / uii

        if i == n
            continue
        end

        ui = u[i, i+1:n] / uii # shape (n-i, )
        xi = x[i+1:n]          # shape (n-i, )
        x[i] -= ui' * xi 
    end

    return BackSubstitutionSolution(x)
end

BackSubstitution = BackSubstitutionProblem
export BackSubstitution
