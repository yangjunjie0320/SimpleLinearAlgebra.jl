abstract type EigenProblemMixin <: ProblemMixin end

module Eigen
    using LinearAlgebra

    import SimpleLinearAlgebra.EigenProblemMixin
    import SimpleLinearAlgebra.SolutionMixin
    
    """
    Base struct for eigenvalue solution that contains eigenvalues and eigenvectors
    """
    struct EigenSolution <: SolutionMixin
        e::AbstractVector # Eigenvalues
        c::AbstractMatrix # Eigenvectors as columns
        is_converged::Bool
    end

    """
    Power Iteration Method - finds the dominant eigenvalue and eigenvector
    """
    struct PowerMethod <: EigenProblemMixin
        a::AbstractMatrix
        max_cycle::Int
        tol::Real
        
        function PowerMethod(a::AbstractMatrix, max_cycle::Int=1000, tol::Real=1e-10)
            @assert size(a, 1) == size(a, 2) "Matrix must be square"
            new(a, max_cycle, tol)
        end
    end

    function kernel(prob::PowerMethod)
        a = prob.a
        n = size(a, 1)
        max_cycle = prob.max_cycle
        tol = prob.tol
        
        # Initialize with random vector
        de = 1.0
        e_pre = 0.0
        e_cur = 0.0

        v_pre = normalize(rand(n))
        v_cur = v_pre
        
        cycle = 1
        is_converged = false
        is_max_cycle = false

        while !is_converged && !is_max_cycle
            is_converged = de < tol
            is_max_cycle = cycle >= max_cycle

            v_cur = a * v_pre
            e_cur = v_pre' * v_cur
            v_cur = normalize(v_cur)

            de = abs(e_cur - e_pre)
            e_pre = e_cur
            v_pre = v_cur
            cycle += 1
        end

        soln = EigenSolution([e_cur], reshape(v_cur, n, 1), is_converged)
        return soln
    end

    struct RayleighQuotient <: EigenProblemMixin
        a::AbstractMatrix
        max_cycle::Int
        tol::Real

        function RayleighQuotient(a::AbstractMatrix, max_cycle::Int=1000, tol::Real=1e-10)
            @assert size(a, 1) == size(a, 2) "Matrix must be square"
            new(a, max_cycle, tol)
        end
    end

    function kernel(prob::RayleighQuotient)
        a = prob.a
        n = size(a, 1)
        max_cycle = prob.max_cycle
        tol = prob.tol
        
        # Initialize with random vector
        de = 1.0
        e_pre = 0.0
        e_cur = 0.0

        v_pre = normalize(rand(n))
        v_cur = v_pre
        
        cycle = 1
        is_converged = false
        is_max_cycle = false

        while !is_converged && !is_max_cycle
            is_converged = de < tol
            is_max_cycle = cycle >= max_cycle

            e_cur = v_pre' * a * v_pre
            v_cur = (a - e_cur * I) \ v_pre
            v_cur = normalize(v_cur)

            de = abs(e_cur - e_pre)
            e_pre = e_cur
            v_pre = v_cur
            cycle += 1
        end

        soln = EigenSolution([e_cur], reshape(v_cur, n, 1), is_converged)
        return soln
    end
    
end

kernel(prob::EigenProblemMixin) = Eigen.kernel(prob)
export Eigen
