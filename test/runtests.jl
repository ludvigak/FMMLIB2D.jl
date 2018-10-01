using FMMLIB2D
using Compat.Test

@testset "Laplace FMM 2D" begin
    include("laplace.jl")
end

@testset "Helmholtz FMM 2D" begin
    include("helmholtz.jl")
end

@testset "Complex FMM 2D" begin
    include("complex.jl")
end
