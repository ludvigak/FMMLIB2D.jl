using FMMLIB2D
using SpecialFunctions
using Compat.LinearAlgebra
using Compat.Random
using Compat.Test

if VERSION < v"0.7"
    srand(0)
else
    Random.seed!(0)
end

function direct_self(source, dipstr)
    pot = zero(source)
    grad = zero(source)
    hess = zero(source)        
    N = length(source)
    for i=1:N
        for j=1:N
            if i == j
                continue
            end
            pot[i] += dipstr[j] / (source[i] - source[j])
            grad[i] += -dipstr[j] / (source[i] - source[j])^2
            hess[i] += 2*dipstr[j] / (source[i] - source[j])^3
        end
    end
    return pot, grad, hess
end

function direct_targ(source, dipstr, target)
    M = length(target)
    N = length(source)
    pot = zero(target)
    grad = zero(target)    
    hess = zero(target)        
    for i=1:M
        for j=1:N
            pot[i] += dipstr[j] / (target[i] - source[j])
            grad[i] += -dipstr[j] / (target[i] - source[j])^2
            hess[i] += 2*dipstr[j] / (target[i] - source[j])^3
        end
    end
    return pot, grad, hess
end    

function runtest()
    iprec = 5
    N = 30
    M = 40
    source = rand(N) + 1im*rand(N)
    dipstr = rand(N) + 1im*rand(N)
    target = rand(M) + 1im*rand(M)
    ifpot = 1
    ifgrad = 1
    ifhess = 1
    ifpottarg = 1
    ifgradtarg = 1
    ifhesstarg = 1

    refpottarg, refgradtarg, refhesstarg = direct_targ(source, dipstr, target)
    refpot, refgrad, refhess = direct_self(source, dipstr)

    @testset "Direct interface" begin
        pot, grad, hess, pottarg, gradtarg, hesstarg = 
            zfmm2dparttarg(iprec, source, dipstr, ifpot, ifgrad, ifhess,
                           target, ifpottarg, ifgradtarg, ifhesstarg)
        
        pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
        grad_relerr = norm(refgrad-grad, Inf) / norm(refgrad, Inf)
        hess_relerr = norm(refhess-hess, Inf) / norm(refhess, Inf)
        
        pottarg_relerr = norm(refpottarg-pottarg, Inf) / norm(refpottarg, Inf)
        gradtarg_relerr = norm(refgradtarg-gradtarg, Inf) / norm(refgradtarg, Inf)
        hesstarg_relerr = norm(refhesstarg-hesstarg, Inf) / norm(refhesstarg, Inf)
        
        @test pot_relerr < 1e-15
        @test grad_relerr < 1e-15
        @test hess_relerr < 1e-15        
        @test pottarg_relerr < 1e-15
        @test gradtarg_relerr < 1e-15
        @test hesstarg_relerr < 1e-15                
    end

    @testset "Keyword interface" begin
        U = zfmm2d(source=source, dipstr=dipstr)
        @test norm(U.pot-refpot, Inf) / norm(refpot, Inf) < 1e-14
        U = zfmm2d(source=source, dipstr=dipstr, target=target)        
        @test norm(U.pottarg-refpottarg, Inf) / norm(refpot, Inf) < 1e-14
    end

    @testset "Tolerance test" begin
        N = 1000
        source = rand(N) + 1im*rand(N)
        dipstr = rand(N) + 1im*rand(N)
        refpot, _, _ = direct_self(source, dipstr)
        for tolexp=-14:-1
            TOL = 0.49*10.0^tolexp
            U = zfmm2d(source=source, tol=TOL, dipstr=dipstr)
            relerr = norm(U.pot-refpot) / norm(refpot)
            @test relerr < TOL
        end
    end
end

runtest()
