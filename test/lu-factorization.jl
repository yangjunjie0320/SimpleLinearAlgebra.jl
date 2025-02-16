@testset "LU factorization version 1" begin
    n = 10
    tol = 1e-10
    a = rand(n, n)

    p = LUFactVersion1(a, tol)
    s = kernel(p)
    l = s.l
    u = s.u

    @test istril(l)
    @test istriu(u)
    @test maximum(abs, a - l * u) < tol
end

@testset "LU factorization version 2" begin
    n = 10
    tol = 1e-10
    a = rand(n, n)

    p = LUFactVersion2(a, tol)
    s = kernel(p)
    l = s.l
    u = s.u

    @test istril(l)
    @test istriu(u)
    @test maximum(abs, a - l * u) < tol
end

@testset "LU factorization" begin
    n = 20
    tol = 1e-10
    a = rand(n, n)

    p = LUFactorizationProblem(a, tol)
    s = kernel(p)
    l = s.l
    u = s.u

    @test istril(l)
    @test istriu(u)
    @test maximum(abs, a - l * u) < tol
end

@testset "singular matrix" begin
    n = 10
    tol = 1e-10
    a = zeros(n, n)
    p = LUFactorizationProblem(a, tol)
    @test_throws LinearAlgebraError kernel(p)
end
