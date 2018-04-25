using FMMLIB2D
using Base.Test

@testset "Laplace FMM 2D" begin
    include("laplace.jl")
end

@testset "Helmholtz FMM 2D" begin
    include("helmholtz.jl")
end
