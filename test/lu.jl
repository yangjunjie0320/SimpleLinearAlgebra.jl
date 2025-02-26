import SimpleLinearAlgebra.LU

@testset "LU factorization version 1" begin
    n = 10
    a = rand(n, n)

    prob = LU.Version1(a)
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

    prob = LU.Version2(a)
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
    n = 10
    a = zeros(n, n)
    prob = SimpleLinearAlgebra.LUFactorization(a)
    @test_throws AssertionError kernel(prob)
end

@testset "Partial pivoting LU factorization" begin
    a = [1e-8 1; 1 1]
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
