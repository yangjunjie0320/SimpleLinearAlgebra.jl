@testset "forward substitution" begin
    n = 10
    tol = 1e-10
    l = LinearAlgebra.tril(rand(n, n))
    b = rand(n)

    p = ForwardSubstitution(l, b, tol)
    s = kernel(p)
    x = s.x
    @test maximum(abs, l * x - b) < tol
end

@testset "not lower triangular" begin
    n = 10
    tol = 1e-10
    l = rand(n, n)
    b = rand(n)
    @test_throws AssertionError ForwardSubstitution(l, b, tol)
end

@testset "singular matrix" begin
    n = 10
    tol = 1e-10
    l = zeros(n, n)
    b = rand(n)
    prob = ForwardSubstitution(l, b, tol)
    @test_throws AssertionError kernel(prob)
end
