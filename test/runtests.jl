using SimpleLinearAlgebra
using Test, LinearAlgebra

@testset "forward substitution" begin
    include("forward.jl")
end

@testset "back substitution" begin
    include("back.jl")
end

@testset "LU factorization" begin
    include("lu.jl")
end

@testset "QR factorization" begin
    include("qr.jl")
end
