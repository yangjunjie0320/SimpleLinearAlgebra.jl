@testset "Partial pivoting LU factorization" begin
    n = 4
    tol = 1e-10
    a = rand(n, n)

    p = PartialPivotingLUFactorization(a, tol)
    s = kernel(p)

    @test istril(s.l)
    @test istriu(s.u)
    @test isperm(s.p)
    
    l = s.l
    u = s.u
    p = zeros(Int, n, n)
    setindex!.(Ref(p), 1, 1:n, s.p)

    @test minimum(a[s.p, :] - l * u) < tol
    @test minimum(p * a - l * u) < tol
end
