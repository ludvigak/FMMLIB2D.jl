__precompile__()
module FMMLIB2D

export lfmm2d, rfmm2d
export lfmm2dparttarg, lfmm2dpartself
export rfmm2dparttarg, rfmm2dpartself

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("FMMLIB2D is not properly installed. Please run Pkg.build(\"FMMLIB2D\").")
end

struct FMMOutput
    pot
    grad
    hess
    pottarg
    gradtarg
    hesstarg
end

## KEYWORD INTERFACES

"""
    lfmm2d(;source, target, charge, dipstr, dipvec, tol,
            ifpot, ifgrad, ifhess,
            ifpottarg, ifgradtarg, ifhesstarg)

    Laplace particle target FMM in R^2 (complex). Keyword interface.
"""
function lfmm2d(;
    source::Array{Float64} = zeros(2, 0),
    target::Array{Float64} = zeros(2, 0),
    charge::Array{Complex128} = zeros(Complex128, 0),
    dipstr::Array{Complex128} = zeros(Complex128, 0),
    dipvec::Array{Float64} = zeros(2, 0),   
    tol::Float64 = 1e-15,
    ifpot::Bool = true,
    ifgrad::Bool = false,
    ifhess::Bool = false,
    ifpottarg::Bool = true,
    ifgradtarg::Bool = false,
    ifhesstarg::Bool = false
                )
    # Parse keywords
    iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
    target, ifpottarg, ifgradtarg, ifhesstarg =
        parse_keywords(tol, source, charge, dipstr, dipvec, ifpot, ifgrad, ifhess,
                       target, ifpottarg, ifgradtarg, ifhesstarg)
    # Call FMM and pack output
    pot, grad, hess, pottarg, gradtarg, hesstarg =        
        lfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot,
                       ifgrad, ifhess, target, ifpottarg, ifgradtarg, ifhesstarg)    
    return FMMOutput(pot, grad, hess, pottarg, gradtarg, hesstarg)
end

"""
    lfmm2d(;source, target, charge, dipstr, dipvec, tol,
            ifpot, ifgrad, ifhess,
            ifpottarg, ifgradtarg, ifhesstarg)

    Laplace particle target FMM in R^2 (real). Keyword interface.
"""
function rfmm2d(;
    source::Array{Float64} = zeros(2, 0),
    target::Array{Float64} = zeros(2, 0),
    charge::Array{Float64} = zeros(Float64, 0),
    dipstr::Array{Float64} = zeros(Float64, 0),
    dipvec::Array{Float64} = zeros(2, 0),   
    tol::Float64 = 1e-15,
    ifpot::Bool = true,
    ifgrad::Bool = false,
    ifhess::Bool = false,
    ifpottarg::Bool = true,
    ifgradtarg::Bool = false,
    ifhesstarg::Bool = false
                )
    # Parse keywords
    iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
    target, ifpottarg, ifgradtarg, ifhesstarg =
        parse_keywords(tol, source, charge, dipstr, dipvec, ifpot, ifgrad, ifhess,
                       target, ifpottarg, ifgradtarg, ifhesstarg)
    # Call FMM and pack output
    pot, grad, hess, pottarg, gradtarg, hesstarg =        
        rfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot,
                       ifgrad, ifhess, target, ifpottarg, ifgradtarg, ifhesstarg)    
    return FMMOutput(pot, grad, hess, pottarg, gradtarg, hesstarg)
end


## DIRECT INTERFACES

"""
    pot, grad, hess, pottarg, gradtarg, hesstarg = 
        lfmm2dparttarg( iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, 
                        ifpot, ifgrad, ifhess, 
                        target, ifpottarg, ifgradtarg, ifhesstarg)

    Laplace particle target FMM in R^2 (complex). Direct library interface.
"""
function lfmm2dparttarg(iprec::Int64,
                        source::Array{Float64},
                        ifcharge::Int64,
                        charge::Array{Complex128},
                        ifdipole::Int64,
                        dipstr::Array{Complex128},
                        dipvec::Array{Float64},
                        ifpot::Int64,
                        ifgrad::Int64,
                        ifhess::Int64,
                        target::Array{Float64},
                        ifpottarg::Int64,
                        ifgradtarg::Int64,
                        ifhesstarg::Int64)
    # Size checks
    nsource = size(source, 2)
    @assert size(source, 1) == 2
    if ifcharge != 0
        @assert length(charge) == nsource
    end
    if ifdipole != 0
        @assert length(dipstr) == nsource
        @assert size(dipvec) == (2, nsource)
    end
    ntarget = size(target, 2)
    @assert size(source, 1) == 2
    # Prepare output structures, only allocate those needed
    ier = 0
    pot = zeros(Complex128, nsource * (ifpot!=0 ? 1 : 0) )
    grad = zeros(Complex128, 2, nsource * (ifgrad!=0 ? 1 : 0) )
    hess = zeros(Complex128, 3, nsource * (ifhess!=0 ? 1 : 0) )
    pottarg = Array{Complex128}(ntarget * (ifpottarg!=0 ? 1 : 0) )
    gradtarg = Array{Complex128}(2, ntarget * (ifgradtarg!=0 ? 1 : 0) )
    hesstarg = Array{Complex128}(3, ntarget * (ifhesstarg!=0 ? 1 : 0) )
    # Library call
    ccall( (:lfmm2dparttarg_, fmmlib2d), Void,
           (Ref{Int64}, # ier
            Ref{Int64}, # iprec
            Ref{Int64}, # nsource
            Ref{Float64}, # source
            Ref{Int64}, # ifcharge
            Ref{Complex128}, #charge
            Ref{Int64}, # ifdipole
            Ref{Complex128}, #dipstr
            Ref{Float64}, #dipvec
            Ref{Int64}, # ifpot
            Ref{Complex128}, # pot
            Ref{Int64}, # ifgrad
            Ref{Complex128}, #grad
            Ref{Int64}, # ifhess
            Ref{Complex128}, #hess
            Ref{Int64}, # ntarget
            Ref{Float64}, # target
            Ref{Int64}, # ifpottarg
            Ref{Complex128}, # pottarg
            Ref{Int64}, # ifgradtarg
            Ref{Complex128}, # gradtarg
            Ref{Int64}, # ifhesstarg
            Ref{Complex128}, # hesstarg
            ),
           ier,iprec,nsource,source,
           ifcharge,charge,ifdipole,dipstr,dipvec,
           ifpot,pot,ifgrad,grad,ifhess,hess,
           ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,ifhesstarg,hesstarg
           )
    # Check and return
    @assert ier == 0    
    return pot, grad, hess, pottarg, gradtarg, hesstarg
end

"""
    pot, grad, hess = 
        lfmm2dpartself( iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, 
                        ifpot, ifgrad, ifhess)

    Laplace particle FMM in R^2 (complex). Direct library interface.
"""
function lfmm2dpartself(iprec::Int64,
                        source::Array{Float64},
                        ifcharge::Int64,
                        charge::Array{Complex128},
                        ifdipole::Int64,
                        dipstr::Array{Complex128},
                        dipvec::Array{Float64},
                        ifpot::Int64,
                        ifgrad::Int64,
                        ifhess::Int64)
    target = zeros(0,0)
    ifpottarg = 0
    ifgradtarg = 0
    ifhesstarg = 0
    pot, grad, hess, _, _, _ = 
        lfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
                       target, ifpottarg, ifgradtarg, ifhesstarg)
    return pot, grad, hess
end

"""
    pot, grad, hess, pottarg, gradtarg, hesstarg = 
        rfmm2dparttarg( iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, 
                        ifpot, ifgrad, ifhess, 
                        target, ifpottarg, ifgradtarg, ifhesstarg)

    Laplace particle target FMM in R^2 (real). Direct library interface.
"""
function rfmm2dparttarg(iprec::Int64,
                        source::Array{Float64},
                        ifcharge::Int64,
                        charge::Array{Float64},
                        ifdipole::Int64,
                        dipstr::Array{Float64},
                        dipvec::Array{Float64},
                        ifpot::Int64,
                        ifgrad::Int64,
                        ifhess::Int64,
                        target::Array{Float64},
                        ifpottarg::Int64,
                        ifgradtarg::Int64,
                        ifhesstarg::Int64)
    # Size checks
    nsource = size(source, 2)
    @assert size(source, 1) == 2
    if ifcharge != 0
        @assert length(charge) == nsource
    end
    if ifdipole != 0
        @assert length(dipstr) == nsource
        @assert size(dipvec) == (2, nsource)
    end
    ntarget = size(target, 2)
    @assert size(source, 1) == 2
    # Prepare output structures, only allocate those needed
    ier = 0
    pot = zeros(Float64, nsource * (ifpot!=0 ? 1 : 0) )
    grad = zeros(Float64, 2, nsource * (ifgrad!=0 ? 1 : 0) )
    hess = zeros(Float64, 3, nsource * (ifhess!=0 ? 1 : 0) )
    pottarg = Array{Float64}(ntarget * (ifpottarg!=0 ? 1 : 0) )
    gradtarg = Array{Float64}(2, ntarget * (ifgradtarg!=0 ? 1 : 0) )
    hesstarg = Array{Float64}(3, ntarget * (ifhesstarg!=0 ? 1 : 0) )
    # Library call
    ccall( (:rfmm2dparttarg_, fmmlib2d), Void,
           (Ref{Int64}, # ier
            Ref{Int64}, # iprec
            Ref{Int64}, # nsource
            Ref{Float64}, # source
            Ref{Int64}, # ifcharge
            Ref{Float64}, #charge
            Ref{Int64}, # ifdipole
            Ref{Float64}, #dipstr
            Ref{Float64}, #dipvec
            Ref{Int64}, # ifpot
            Ref{Float64}, # pot
            Ref{Int64}, # ifgrad
            Ref{Float64}, #grad
            Ref{Int64}, # ifhess
            Ref{Float64}, #hess
            Ref{Int64}, # ntarget
            Ref{Float64}, # target
            Ref{Int64}, # ifpottarg
            Ref{Float64}, # pottarg
            Ref{Int64}, # ifgradtarg
            Ref{Float64}, # gradtarg
            Ref{Int64}, # ifhesstarg
            Ref{Float64}, # hesstarg
            ),
           ier,iprec,nsource,source,
           ifcharge,charge,ifdipole,dipstr,dipvec,
           ifpot,pot,ifgrad,grad,ifhess,hess,
           ntarget,target,ifpottarg,pottarg,ifgradtarg,gradtarg,ifhesstarg,hesstarg
           )
    # Check and return
    @assert ier == 0    
    return pot, grad, hess, pottarg, gradtarg, hesstarg
end

"""
    pot, grad, hess = 
        rfmm2dpartself( iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, 
                        ifpot, ifgrad, ifhess)

    Laplace particle FMM in R^2 (real). Direct library interface.
"""
function rfmm2dpartself(iprec::Int64,
                        source::Array{Float64},
                        ifcharge::Int64,
                        charge::Array{Float64},
                        ifdipole::Int64,
                        dipstr::Array{Float64},
                        dipvec::Array{Float64},
                        ifpot::Int64,
                        ifgrad::Int64,
                        ifhess::Int64)
    target = zeros(0,0)
    ifpottarg = 0
    ifgradtarg = 0
    ifhesstarg = 0
    pot, grad, hess, _, _, _ = 
        rfmm2dparttarg(iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad, ifhess,
                       target, ifpottarg, ifgradtarg, ifhesstarg)
    return pot, grad, hess
end

## HELPER FUNCTIONS

function parse_keywords(tol, source, charge, dipstr, dipvec, ifpot, ifgrad, ifhess,
                        target, ifpottarg, ifgradtarg, ifhesstarg)
    # Set correct integer flags
    if length(target)==0
        ifpottarg = false
        ifgradtarg = false
        ifhesstarg = false
    end
    ifcharge = length(charge)>0 ? 1 : 0
    ifdipole = length(dipstr)>0 ? 1 : 0
    ifpot = ifpot ? 1 : 0
    ifgrad = ifgrad ? 1 : 0
    ifhess = ifhess ? 1 : 0
    ifpottarg = ifpottarg ? 1 : 0
    ifgradtarg = ifgradtarg ? 1 : 0
    ifhesstarg = ifhesstarg ? 1 : 0    
    # Select the right precision flag
    iprec = tol2iprec(tol)
    return iprec, source, ifcharge, charge, ifdipole, dipstr, dipvec, ifpot, ifgrad,
    ifhess, target, ifpottarg, ifgradtarg, ifhesstarg
end

function tol2iprec(tol)
    if tol < 0.5e-14
        iprec = 5
    elseif tol < 0.5e-12
        iprec = 4
    elseif tol < 0.5e-9
        iprec = 3
    elseif tol < 0.5e-6
        iprec = 2
    elseif tol < 0.5e-3
        iprec = 1
    elseif tol < 0.5e-2
        iprec =
            0
    elseif tol < 0.5e-1
        iprec = -1
    else
        iprec = -2
    end
end

end # module
