using FMMLIB2D
using Base.Test

srand(0)

# TODO: Include dipoles
# TODO: Compute gradient and hessian

function direct_self(source, charge)
    refpot = zeros(charge)
    N = size(source, 2)
    for i=1:N
        for j=1:N
            if i == j
                continue
            end
            refpot[i] += charge[j]*log(norm(source[:,i]-source[:,j]))
        end
    end
    return refpot
end
function direct_targ(source, charge, target)
    M = size(target, 2)
    N = size(source, 2)
    refpottarg = zeros(typeof(charge[1]), M)    
    for i=1:M
        for j=1:N
            refpottarg[i] += charge[j]*log(norm(target[:,i]-source[:,j]))
        end
    end
    return refpottarg
end    

@testset "Laplace FMM 2D (complex)" begin
    iprec = 5
    N = 40
    M = 30
    source = rand(2, N)
    charge = rand(N) + 1im*rand(N)
    dipstr = complex(rand(N))
    dipvec = rand(2, N)
    target = rand(2, M)
    ifcharge = 1
    ifdipole = 0
    ifpot = 1
    ifgrad = 1
    ifhess = 1
    ifpottarg = 1
    ifgradtarg = 1
    ifhesstarg = 1
    refpottarg = direct_targ(source, charge, target)
    refpot = direct_self(source, charge)
    @testset "Direct interface: lfmm2dpartself" begin
        pot, grad, hess = 
            lfmm2dpartself(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess)
        pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
        @test pot_relerr < 1e-15
    end
    @testset "Direct interface: lfmm2dparttarg" begin
        pot, grad, hess, pottarg, gradtarg, hesstarg = 
            lfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
                           target, ifpottarg, ifgradtarg, ifhesstarg)
        pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
        pottarg_relerr = norm(refpottarg-pottarg, Inf) / norm(refpottarg, Inf)    
        @test pot_relerr < 1e-15
        @test pottarg_relerr < 1e-15
    end
    @testset "Keyword interface" begin
        U = lfmm2d(source=source, charge=charge, target=target)
        @test norm(U.pot-refpot, Inf) / norm(refpot, Inf) < 1e-14
        @test norm(U.pottarg-refpottarg, Inf) / norm(refpot, Inf) < 1e-14
    end    
end

@testset "Laplace FMM 2D (real)" begin
    iprec = 5
    N = 40
    M = 30
    source = rand(2, N)
    charge = rand(N)
    dipstr = rand(N)
    dipvec = rand(2, N)
    target = rand(2, M)
    ifcharge = 1
    ifdipole = 0
    ifpot = 1
    ifgrad = 1
    ifhess = 1
    ifpottarg = 1
    ifgradtarg = 1
    ifhesstarg = 1
    refpottarg = direct_targ(source, charge, target)
    refpot = direct_self(source, charge)
    @testset "Direct interface: rfmm2dpartself" begin
        pot, grad, hess = 
            rfmm2dpartself(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess)
        pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
        @test pot_relerr < 1e-15
    end    
    @testset "Direct interface: rfmm2dparttarg" begin
        pot, grad, hess, pottarg, gradtarg, hesstarg = 
            rfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
                           target, ifpottarg, ifgradtarg, ifhesstarg)
        pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
        pottarg_relerr = norm(refpottarg-pottarg, Inf) / norm(refpottarg, Inf)    
        @test pot_relerr < 1e-15
        @test pottarg_relerr < 1e-15
    end
    @testset "Keyword interface" begin
        U = rfmm2d(source=source, charge=charge, target=target)
        @test norm(U.pot-refpot, Inf) / norm(refpot, Inf) < 1e-14
        @test norm(U.pottarg-refpottarg, Inf) / norm(refpot, Inf) < 1e-14
    end    
end
