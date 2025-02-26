# QR factorization
# A = QR
#
# where Q is orthogonal and R is upper triangular

@testset "QR factorization" begin
    n = 10
    tol = 1e-10
    a = rand(n, n)

    p = QRFactorization(a)
    s = kernel(p)

    q = s.q
    r = s.r

    @test istriu(r)
    @test q * q' ≈ Matrix{Float64}(I, n, n)
    @test minimum(a - q * r) < tol
end
