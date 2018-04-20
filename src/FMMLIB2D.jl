module FMMLIB2D

export lfmm2dpartself, lfmm2dparttarg

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("FMMLIB2D is not properly installed. Please run Pkg.build(\"FMMLIB2D\").")
end

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

end # module
