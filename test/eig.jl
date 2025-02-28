import SimpleLinearAlgebra.Eigen

@testset "Power Iteration Method" begin
    # Create a matrix with a known dominant eigenvalue
    n = 8
    a = randn(n, n)
    # Make the first eigenvalue significantly larger
    a = a * a'  # Create a symmetric positive-definite matrix

    prob = Eigen.PowerMethod(a, 1000, 1e-10)
    tol = sqrt(prob.tol) * n
    
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

    prob = Eigen.RayleighQuotient(a, 1000, 1e-10)
    tol = sqrt(prob.tol) * n
    
    soln = kernel(prob)
    e = soln.e[1]
    c = soln.c[:, 1]

    @test isapprox(a * c, e * c, atol=tol)
end

@testset "Singular Value Decomposition" begin
    n = 4
    a = randn(n, n)
    a = a * a'  # Create a symmetric positive-definite matrix

    prob = SingularValueDecomposition(a, 1e-8)
    tol = prob.tol

    soln = kernel(prob)
    u = soln.u
    s = soln.s
    v = soln.v

    @test isapprox(a, u * s * v', atol=tol)
end