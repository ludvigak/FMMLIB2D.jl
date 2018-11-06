# FMMLIB2D 
[![Build Status](https://travis-ci.org/ludvigak/FMMLIB2D.jl.svg?branch=master)](https://travis-ci.org/ludvigak/FMMLIB2D.jl)
[![Coverage Status](https://coveralls.io/repos/github/ludvigak/FMMLIB2D.jl/badge.svg?branch=master)](https://coveralls.io/github/ludvigak/FMMLIB2D.jl?branch=master)

This is a Julia interface to the Fast Multipole Method (FMM) library
[FMMLIB2D](https://github.com/zgimbutas/fmmlib2d) by Leslie Greengard and Zydrunas
Gimbutas.

Documentation for the library can be found in the [FMMLIB2D User's Guide](https://github.com/ludvigak/fmmlib2d/blob/master/doc/fmm2dpart_manual.pdf).

This package currently provides interfaces to the FMM's for Laplace (real and complex), Helmholtz, and complex sums: 
`rfmm2dpartself`, `rfmm2dparttarg`, `lfmm2dpartself`, `lfmm2dparttarg`, `hfmm2dparttarg`, `zfmm2dparttarg`

The most convenient way of calling them is through the Julia interfaces with keyword arguments, e.g. 
```julia
x = rand(2, 10)
y = rand(2, 20)
q = rand(10) + 1im*rand(10)
U = lfmm2d(source=x, charge=q, target=y, ifgradtarg=true, tol=1e-9)
```

### Real Laplace FMM:
```julia
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
```

Output format:
```julia
U.pot      (Nsrc)
U.grad     (2,Nsrc)
U.hess     (3,Nsrc)
U.pottarg  (Ntrg)
U.gradtarg (2,Ntrg)
U.hesstarg (3,Ntrg)
```


### Complex Laplace FMM:
```julia
U = lfmm2d(source::Array{Float64} = ...,
           target::Array{Float64} = ...,
           charge::Array{ComplexF64} = ...,
           dipstr::Array{ComplexF64} = ...,
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

Output format:
```julia
U.pot      (Nsrc)
U.grad     (2,Nsrc)
U.hess     (3,Nsrc)
U.pottarg  (Ntrg)
U.gradtarg (2,Ntrg)
U.hesstarg (3,Ntrg)
```

### Helmholtz FMM:
```julia
U = hfmm2d(source::Array{Float64} = ...,
           target::Array{Float64} = ...,
           charge::Array{ComplexF64} = ...,
           dipstr::Array{ComplexF64} = ...,
           dipvec::Array{Float64} = ...,
           tol::Float64 = ...,
           zk::ComplexF64 = ...,
           ifpot::Bool = ...,
           ifgrad::Bool = ...,
           ifhess::Bool = ...,
           ifpottarg::Bool = ...,
           ifgradtarg::Bool = ...,
           ifhesstarg::Bool = ...,
           )
```

Output format:
```julia
U.pot      (Nsrc)
U.grad     (2,Nsrc)
U.hess     (3,Nsrc)
U.pottarg  (Ntrg)
U.gradtarg (2,Ntrg)
U.hesstarg (3,Ntrg)
```

### Complex FMM:

```julia
U = zfmm2d(source::Array{ComplexF64} = ...,
           target::Array{ComplexF64} = ...,
           dipstr::Array{ComplexF64} = ...,
           tol::Float64 = 1e-15,
           ifpot::Bool = true,
           ifgrad::Bool = false,
           ifhess::Bool = false,
           ifpottarg::Bool = true,
           ifgradtarg::Bool = false,
           ifhesstarg::Bool = false
           )
```

Output format:
```julia
U.pot      (Nsrc)
U.grad     (Nsrc)
U.hess     (Nsrc)
U.pottarg  (Ntrg)
U.gradtarg (Ntrg)
U.hesstarg (Ntrg)
```
