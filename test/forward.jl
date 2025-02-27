using LinearAlgebra

@testset "forward substitution" begin
    n = 10
    l = tril(rand(n, n))
    b = rand(n)

    prob = ForwardSubstitution(l, b)
    tol = prob.tol

    soln = kernel(prob)
    x = soln.x

    @test isapprox(l * x, b, atol=tol)
end

@testset "not lower triangular" begin
    n = 10
    l = rand(n, n)
    b = rand(n)
    @test_throws AssertionError ForwardSubstitution(l, b)
end

@testset "singular matrix" begin
    n = 10
    l = zeros(n, n)
    b = rand(n)
    prob = ForwardSubstitution(l, b)
    @test_throws AssertionError kernel(prob)
end
