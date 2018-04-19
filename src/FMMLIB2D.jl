module FMMLIB2D

export lfmm2dpartself

const depsfile = joinpath(dirname(@__DIR__), "deps", "deps.jl")
if isfile(depsfile)
    include(depsfile)
else
    error("FMMLIB2D is not properly installed. Please run Pkg.build(\"FMMLIB2D\").")
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
    nsource = size(source, 2)
    @assert size(source, 1) == 2
    if ifcharge != 0
        @assert length(charge) == nsource
    end
    if ifdipole != 0
        @assert length(dipstr) == nsource
        @assert size(dipvec) == (2, nsource)
    end    
    ier = zeros(Int64, 1)
    pot = zeros(Complex128, nsource)
    grad= zeros(Complex128, 2, nsource)
    hess = zeros(Complex128, 3, nsource)
    pottarg = Array{Complex128}(1)
    gradtarg = Array{Complex128}(2,1)
    hesstarg = Array{Complex128}(3,1)    
    ccall( (:lfmm2dpartself_, fmmlib2d), Void,
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
            Ref{Complex128}, #pot
            Ref{Int64}, # ifgrad
            Ref{Complex128}, #grad
            Ref{Int64}, # ifhess
            Ref{Complex128} #hess
            ),
           ier,iprec,nsource,source,
           ifcharge,charge,ifdipole,dipstr,dipvec,
           ifpot,pot,ifgrad,grad,ifhess,hess           
           )
    ier = ier[1]
    @assert ier == 0    
    return pot, grad, hess
end

end # module
