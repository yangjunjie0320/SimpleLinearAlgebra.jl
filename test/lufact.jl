@testset "LU factorization" begin
    n = 10
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
