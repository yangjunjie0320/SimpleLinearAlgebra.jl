import SimpleLinearAlgebra.Eigen

@testset "Power Iteration Method" begin
    # Create a matrix with a known dominant eigenvalue
    n = 4
    a = randn(n, n)
    # Make the first eigenvalue significantly larger
    a = a * a'  # Create a symmetric positive-definite matrix
    eig_true = eigen(a)

    println(eig_true.values)
    
    prob = Eigen.PowerIteration(a, 1000, 1e-8)
    tol = sqrt(prob.tol)
    
    soln = kernel(prob)
    e = soln.e[1]
    c = soln.c[:, 1]
    
    # Check eigenvector property: Ax = λx
    if soln.is_converged
        @test isapprox(a * c, e * c, atol=tol)
    end
end
