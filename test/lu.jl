import SimpleLinearAlgebra.LU

@testset "LU factorization version 1" begin
    n = 10
    a = rand(n, n)

    prob = LU.GaussianEliminationV1(a)
    tol = prob.tol

    soln = kernel(prob)
    l = soln.l
    u = soln.u
    
    @test istril(l)
    @test istriu(u)
    @test isapprox(a, l * u, atol=tol)
end

@testset "LU factorization version 2" begin
    n = 10
    a = rand(n, n)

    prob = LU.GaussianEliminationV2(a)
    tol = prob.tol

    soln = kernel(prob)
    l = soln.l
    u = soln.u

    @test istril(l)
    @test istriu(u)
    @test isapprox(a, l * u, atol=tol)
end

@testset "LU factorization" begin
    n = 20
    a = rand(n, n)

    prob = LUFactorization(a)
    tol = prob.tol

    soln = kernel(prob)
    l = soln.l
    u = soln.u

    @test istril(l)
    @test istriu(u)
    @test isapprox(a, l * u, atol=tol)
end

@testset "singular matrix" begin
    a = [0.0 1.0; 1.0 1.0]
    n = size(a, 1)

    prob = LUFactorization(a)
    @test_throws AssertionError kernel(prob)
end

@testset "Partial pivoting LU factorization" begin
    a = [0.0 1.0; 1.0 1.0]
    n = size(a, 1)

    prob = LU.PartialPivoting(a)
    tol = prob.tol

    soln = kernel(prob)
    @test isperm(soln.p)

    l = soln.l
    u = soln.u
    @test istril(l)
    @test istriu(u)
    
    p = zeros(Int, n, n)
    setindex!.(Ref(p), 1, 1:n, soln.p)

    pa = a[soln.p, :]
    @test isapprox(pa, l * u, atol=tol)  # Verify PA = LU
    @test isapprox(p * a, l * u, atol=tol)
end

@testset "Cholesky factorization" begin
    n = 10
    q, r = rand(n, n) |> qr
    a = q * Diagonal(rand(n)) * q'

    prob = CholeskyFactorization(a)
    tol = prob.tol

    soln = kernel(prob)
    l = soln.l
    
    @test istril(l)
    @test isapprox(a, l * l', atol=tol) 
end
