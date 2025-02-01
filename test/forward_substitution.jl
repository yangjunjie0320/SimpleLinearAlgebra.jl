@testset "forward substitution" begin
    n = 10
    tol = 1e-10
    u = LinearAlgebra.tril(rand(n, n))
    b = rand(n)

    p = ForwardSubstitution(u, b, tol)
    s = kernel(p)
    x = s.x
    @test maximum(abs, u * x - b) < tol
end
