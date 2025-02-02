using SimpleLinearAlgebra
using Test, LinearAlgebra

@testset "forward substitution" begin
    include("forward_substitution.jl")
end

@testset "backward substitution" begin
    include("backward_substitution.jl")
end

@testset "LU factorization" begin
    include("lufact.jl")
end
