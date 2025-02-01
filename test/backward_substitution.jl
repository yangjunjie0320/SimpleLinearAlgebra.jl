@testset "backward substitution" begin
    n = 10
    tol = 1e-10
    u = LinearAlgebra.triu(rand(n, n))
    b = rand(n)

    p = BackSubstitution(u, b, tol)
    s = kernel(p)
    x = s.x
    @test maximum(abs, u * x - b) < tol
end
