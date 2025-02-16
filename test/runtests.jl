using SimpleLinearAlgebra
import SimpleLinearAlgebra: LinearAlgebraError

using Test, LinearAlgebra

@testset "forward substitution" begin
    include("forward-substitution.jl")
end

@testset "back substitution" begin
    include("back-substitution.jl")
end

@testset "LU factorization" begin
    include("lu-factorization.jl")
end

@testset "Partial pivoting LU factorization" begin
    include("partial-pivoting-lu-factorization.jl")
end
