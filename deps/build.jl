using BinDeps
using Compat.Libdl

@BinDeps.setup

fmmlib2d = library_dependency("fmmlib2d")

provides(Sources, URI("https://github.com/ludvigak/fmmlib2d/archive/master.zip"),
         fmmlib2d, unpacked_dir = "fmmlib2d-master")

fmmsrcdir = joinpath(BinDeps.srcdir(fmmlib2d), "fmmlib2d-master/src")
libname = "fmmlib2d." * Libdl.dlext
libfile = joinpath(BinDeps.libdir(fmmlib2d),libname)
buildfile = joinpath(fmmsrcdir, libname)

fsrcs = "hfmm2dpart.f hfmm2drouts.f d2tstrcr_omp.f d2mtreeplot.f h2dterms.f helmrouts2d.f cdjseval2d.f hank103.f prini.f cfmm2dpart.f zfmm2dpart.f rfmm2dpart.f lfmm2dpart.f lfmm2drouts.f l2dterms.f laprouts2d.f"

fsrcs = split(fsrcs, " ")
buildcmd = `gfortran -O2 -fPIC -shared $fsrcs -o $libname`

provides(BuildProcess,
         (@build_steps begin
          GetSources(fmmlib2d)
          @build_steps begin
          ChangeDirectory(fmmsrcdir)
          FileRule(libfile, @build_steps begin
                   buildcmd
                   CreateDirectory(libdir(fmmlib2d))
                   `cp $libname $libfile`
                   end)
          end
          end),
         fmmlib2d)

@BinDeps.install Dict(:fmmlib2d => :fmmlib2d)
