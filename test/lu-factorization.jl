@testset "LU factorization version 1" begin
    n = 10
    tol = 1e-10
    a = rand(n, n)

    p = LUFactorizationV1(a, tol)
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

    p = LUFactorizationV2(a, tol)
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

    p = LUFactorization(a, tol)
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
    p = LUFactorization(a, tol)
    @test_throws AssertionError kernel(p)
end
