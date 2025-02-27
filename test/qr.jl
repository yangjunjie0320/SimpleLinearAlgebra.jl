import SimpleLinearAlgebra.QR

@testset "QR factorization with Householder reflections" begin
    n = 10
    tol = 1e-10
    a = rand(n, n)

    prob = QR.HouseholderReflection(a)
    tol = prob.tol

    soln = kernel(prob)
    q = soln.q
    r = soln.r

    @test istriu(r)
    @test isapprox(q' * q, I, atol=tol)
    @test isapprox(a, q * r, atol=tol)
end

@testset "QR factorization with given rotation" begin
    n = 10
    a = rand(n, n)

    prob = QR.GivenRotation(a)
    tol = prob.tol

    soln = kernel(prob)
    q = soln.q
    r = soln.r

    @test istriu(r)
    @test isapprox(q' * q, I, atol=tol)
    @test isapprox(a, q * r, atol=tol)
end

@testset "QR factorization with Gram-Schmidt" begin
    n = 20
    a = rand(n, n)

    prob = QR.GramSchmidt(a)
    tol = prob.tol

    soln = kernel(prob)
    q = soln.q
    r = soln.r

    @test istriu(r)
    @test isapprox(q' * q, I, atol=tol)
    @test isapprox(a, q * r, atol=tol)
end

@testset "QR factorization with modified Gram-Schmidt" begin
    n = 20
    a = rand(n, n)

    prob = QR.ModifiedGramSchmidt(a)
    tol = prob.tol

    soln = kernel(prob)
    q = soln.q
    r = soln.r

    @test istriu(r)
    @test isapprox(q' * q, I, atol=tol)
    @test isapprox(a, q * r, atol=tol)
end

@testset "Symmetric QR factorization" begin
    n = 4
    a = rand(n, n)
    a = a * a'

    prob = SymmetricQRFactorization(a, 1e-6)
    tol = prob.tol

    soln = kernel(prob)
    q = soln.q
    r = soln.r
    
    @test isapprox(r, triu(r, -2), atol=tol)
    @test isapprox(r, tril(r,  2), atol=tol)
    @test isapprox(q' * q, I, atol=tol)
    @test isapprox(a, q * r * q', atol=tol)
end