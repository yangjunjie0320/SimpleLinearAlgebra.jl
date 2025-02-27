import SimpleLinearAlgebra.Eigen

@testset "Power Iteration Method" begin
    # Create a matrix with a known dominant eigenvalue
    n = 8
    a = randn(n, n)
    # Make the first eigenvalue significantly larger
    a = a * a'  # Create a symmetric positive-definite matrix

    prob = Eigen.PowerMethod(a, 1000, 1e-8)
    tol = sqrt(prob.tol * n)
    
    soln = kernel(prob)
    e = soln.e[1]
    c = soln.c[:, 1]

    @test isapprox(a * c, e * c, atol=tol)
end

@testset "Rayleigh Quotient Method" begin
    # Create a matrix with a known dominant eigenvalue
    n = 8
    a = randn(n, n)
    # Make the first eigenvalue significantly larger
    a = a * a'  # Create a symmetric positive-definite matrix

    prob = Eigen.RayleighQuotient(a, 1000, 1e-8)
    tol = sqrt(prob.tol * n)
    
    soln = kernel(prob)
    e = soln.e[1]
    c = soln.c[:, 1]

    @test isapprox(a * c, e * c, atol=tol)
end
