# FMMLIB2D 
[![Build Status](https://travis-ci.org/ludvigak/FMMLIB2D.jl.svg?branch=master)](https://travis-ci.org/ludvigak/FMMLIB2D.jl)
[![Coverage Status](https://coveralls.io/repos/github/ludvigak/FMMLIB2D.jl/badge.svg?branch=master)](https://coveralls.io/github/ludvigak/FMMLIB2D.jl?branch=master)
[![FMMLIB2D](http://pkg.julialang.org/badges/FMMLIB2D_0.6.svg)](http://pkg.julialang.org/detail/FMMLIB2D)

This is a Julia interface to the Fast Multipole Method (FMM) library
[FMMLIB2D](https://github.com/zgimbutas/fmmlib2d) by Leslie Greengard and Zydrunas
Gimbutas.

Documentation for the library can be found in the [FMMLIB2D User's Guide](https://github.com/ludvigak/fmmlib2d/blob/master/doc/fmm2dpart_manual.pdf).

This package currently provides interfaces to the FMM's for Laplace (real and complex) and Helmholtz: 
`rfmm2dpartself`, `rfmm2dparttarg`, `lfmm2dpartself`, `lfmm2dparttarg`, `hfmm2dparttarg`

The most convenient way of calling them are through the Julia interfaces with keyword arguments:

Real Laplace FMM:
```
U = rfmm2d(source::Array{Float64} = ...,
           target::Array{Float64} = ...,
           charge::Array{Float64} = ...,
           dipstr::Array{Float64} = ...,
           dipvec::Array{Float64} = ...,
           tol::Float64 = ...,
           ifpot::Bool = ...,
           ifgrad::Bool = ...,
           ifhess::Bool = ...,
           ifpottarg::Bool = ...,
           ifgradtarg::Bool = ...,
           ifhesstarg::Bool = ...,
           )
Example:
U = fmm2d(source=x, charge=q, target=y, ifpottarg=true, tol=1e-9)
```

Complex Laplace FMM:
```
U = lfmm2d(source::Array{Float64} = ...,
           target::Array{Float64} = ...,
           charge::Array{Complex128} = ...,
           dipstr::Array{Complex128} = ...,
           dipvec::Array{Float64} = ...,
           tol::Float64 = ...,
           ifpot::Bool = ...,
           ifgrad::Bool = ...,
           ifhess::Bool = ...,
           ifpottarg::Bool = ...,
           ifgradtarg::Bool = ...,
           ifhesstarg::Bool = ...,
           )
```

Helmholtz FMM:
```
U = hfmm2d(source::Array{Float64} = ...,
           target::Array{Float64} = ...,
           charge::Array{Complex128} = ...,
           dipstr::Array{Complex128} = ...,
           dipvec::Array{Float64} = ...,
           tol::Float64 = ...,
           zk::Complex128 = ...,
           ifpot::Bool = ...,
           ifgrad::Bool = ...,
           ifhess::Bool = ...,
           ifpottarg::Bool = ...,
           ifgradtarg::Bool = ...,
           ifhesstarg::Bool = ...,
           )
```

Output format:
```
U.pot      (Nsrc)
U.grad     (2,Nsrc)
U.hess     (3,Nsrc)
U.pottarg  (Ntrg)
U.gradtarg (2,Ntrg)
U.hesstarg (3,Ntrg)
```
