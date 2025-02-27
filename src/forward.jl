struct ForwardSubstitutionProblem <: ProblemMixin
    l::AbstractMatrix
    b::AbstractVector
    tol::Real

    function ForwardSubstitutionProblem(l, b, tol=1e-10)
        prob = new(l, b, tol)
        @assert istril(l) "matrix must be lower triangular"
        return prob
    end
end

struct ForwardSubstitutionSolution <: SolutionMixin
    x::AbstractVector
end

function kernel(prob::ForwardSubstitutionProblem)
    l = prob.l
    @assert istril(l) "matrix must be lower triangular"

    b = copy(prob.b)
    n = length(b)
    tol = prob.tol

    x = zeros(n)

    for i in 1:n
        lii = l[i, i]
        @assert abs(lii) > tol "matrix is singular"
        x[i] += b[i] / lii

        if i == 1
            continue
        end
        
        li = l[i, 1:i-1] / lii # shape (i-1, )
        xi = x[1:i-1]          # shape (i-1, )
        x[i] -= li' * xi 
    end
    return ForwardSubstitutionSolution(x)
end

ForwardSubstitution = ForwardSubstitutionProblem
export ForwardSubstitution
