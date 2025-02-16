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

# test for me error message
@testset "not upper triangular" begin
    n = 10
    tol = 1e-10
    u = rand(n, n)
    b = rand(n)
    @test_throws AssertionError BackSubstitution(u, b, tol)
end

@testset "singular matrix" begin
    n = 10
    tol = 1e-10
    u = zeros(n, n)
    b = rand(n)
    prob = BackSubstitution(u, b, tol)
    @test_throws AssertionError kernel(prob)
end
