using LinearAlgebra

@testset "backward substitution" begin
    n = 10
    u = triu(rand(n, n))
    b = rand(n)

    prob = BackSubstitution(u, b)
    tol = prob.tol

    soln = kernel(prob)
    x = soln.x

    @test isapprox(u * x, b, atol=tol)
end

# test for me error message
@testset "not upper triangular" begin
    n = 10
    u = rand(n, n)
    b = rand(n)
    @test_throws AssertionError BackSubstitution(u, b)
end

@testset "singular matrix" begin
    n = 10
    u = zeros(n, n)
    b = rand(n)
    prob = BackSubstitution(u, b)
    @test_throws AssertionError kernel(prob)
end
