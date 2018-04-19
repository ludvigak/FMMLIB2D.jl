using FMMLIB2D
using Base.Test

srand(0)

iprec = 5

N = 5
source = rand(2, N)
charge = rand(N) + 1im*rand(N)
dipstr = complex(rand(N))
dipvec = rand(2, N)

ifcharge = 1
ifdipole = 0
ifpot = 1
ifgrad = 1
ifhess = 1


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