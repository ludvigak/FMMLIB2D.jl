using FMMLIB2D
using Base.Test

srand(0)

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

@testset "Laplace self" begin
    pot, grad, hess = 
        lfmm2dpartself(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess)
    refpot = zeros(pot)
    for i=1:N
        for j=1:N
            if i == j
                continue
            end
            refpot[i] += charge[j]*log(norm(source[:,i]-source[:,j]))
        end
    end
    pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
    @test pot_relerr < 1e-15
end

@testset "Laplace targ" begin
    pot, grad, hess, pottarg, gradtarg, hesstarg = 
        lfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
                       target, ifpottarg, ifgradtarg, ifhesstarg)
    refpot = zeros(pot)    
    for i=1:N
        for j=1:N
            if i == j
                continue
            end
            refpot[i] += charge[j]*log(norm(source[:,i]-source[:,j]))
        end
    end
    pot_relerr = norm(refpot-pot, Inf) / norm(refpot, Inf)
    @test pot_relerr < 1e-15
    refpottarg = zeros(pottarg)
    for i=1:M
        for j=1:N
            refpottarg[i] += charge[j]*log(norm(target[:,i]-source[:,j]))
        end
    end
    pottarg_relerr = norm(refpottarg-pottarg, Inf) / norm(refpottarg, Inf)
    @test pottarg_relerr < 1e-15
end
