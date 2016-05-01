# ApproxMinimumDegree.jl

* [Description](#description)
* [Usage demo](#usage-demo)
* [Perf demo](#performance-demo)
* [Future](#future)
* [Refs](#references)

## Description
ApproxMinimumDegree.jl is a native-Julia implementation of the approximate minimum external
degree algorithm. This is a prototype; prefer [AMD](http://faculty.cse.tamu.edu/davis/suitesparse.html)
([AMD.jl](https://github.com/JuliaOptimizers/AMD.jl)) or [HSL-MC47](http://www.hsl.rl.ac.uk/catalogue/mc47.html)
for anything beyond experimentation.

## Usage demo

```julia
julia> import ApproxMinimumDegree: amd
julia> N = 9;
julia> A = sparse(SymTridiagonal(2*ones(Int, N), -1*ones(Int, N-1)));
julia> rp = randperm(N);
julia> rpA = A[rp, rp];
julia> pamd = amd(rpA);

julia> full(rpA)
9x9 Array{Int64,2}:
  2   0   0   0  -1   0  -1   0   0
  0   2  -1   0   0   0   0  -1   0
  0  -1   2  -1   0   0   0   0   0
  0   0  -1   2   0  -1   0   0   0
 -1   0   0   0   2   0   0   0  -1
  0   0   0  -1   0   2   0   0   0
 -1   0   0   0   0   0   2  -1   0
  0  -1   0   0   0   0  -1   2   0
  0   0   0   0  -1   0   0   0   2
  
julia> full(rpA[pamd, pamd])
9x9 Array{Int64,2}:
  2  -1   0   0   0   0   0   0   0
 -1   2  -1   0   0   0   0   0   0
  0  -1   2  -1   0   0   0   0   0
  0   0  -1   2  -1   0   0   0   0
  0   0   0  -1   2  -1   0   0   0
  0   0   0   0  -1   2  -1   0   0
  0   0   0   0   0  -1   2  -1   0
  0   0   0   0   0   0  -1   2  -1
  0   0   0   0   0   0   0  -1   2
```
`ApproxMinDegree.amd(...)` presently returns a tuple `(perm, iperm)` where `perm` is
an assembly-tree postordering and `iperm` is the corresponding inverse permutation.

## Performance demo

Preordering quality:

```julia
julia> using AMD
julia> using ApproxMinimumDegree
julia> begin

       "Five-point-stencil finite-difference approximation of the bivariate Laplacian."
       function laplacian(N)
          thediag = speye(N, N)
          offdiag = sparse(2:N, 1:(N-1), 1, N, N)
          tridiag = 2*thediag - offdiag - offdiag'
          kron(tridiag, thediag) + kron(thediag, tridiag)
       end

       A = laplacian(30); # matrix order N^2 = 900
       pnative = ApproxMinimumDegree.amd(A);
       pwrapped = AMD.amd(A);

       end;

julia> cholfact(A, perm = 1:30^2)
Base.SparseMatrix.CHOLMOD.Factor{Float64}
type:          LLt
method: simplicial
maxnnz:      27029
nnz:         27029

julia> cholfact(A, perm = pwrapped)
Base.SparseMatrix.CHOLMOD.Factor{Float64}
type:          LLt
method: simplicial
maxnnz:      10231
nnz:         10231

julia> cholfact(A, perm = pnative)
Base.SparseMatrix.CHOLMOD.Factor{Float64}
type:          LLt
method: simplicial
maxnnz:      10082
nnz:         10082
```

Runtime:

```julia
julia> # julia -O --check-bounds=no
julia> versioninfo()
Julia Version 0.4.3
Commit a2f713d* (2016-01-12 21:37 UTC)
Platform Info:
  System: Darwin (x86_64-apple-darwin14.5.0)
  CPU: Intel(R) Core(TM) i7-3520M CPU @ 2.90GHz
  WORD_SIZE: 64
  BLAS: libopenblas (DYNAMIC_ARCH NO_AFFINITY Sandybridge)
  LAPACK: libopenblas
  LIBM: libopenlibm
  LLVM: libLLVM-3.3
julia> using AMD
julia> using Benchmarks
julia> using ApproxMinimumDegree
julia> begin

       "Pretty-print benchmark mean evaluation time and 95% CIs."
       function prettytimes(bench)
           stats = Benchmarks.SummaryStatistics(bench)
           timecenter = stats.elapsed_time_center
           timelower = get(stats.elapsed_time_lower)
           timeupper = get(stats.elapsed_time_upper)
           # based on Benchmarks.pretty_time_string
           tscale, tunits =
               timecenter < 10^4 ? (10^0, "ns") :
               timecenter < 10^7 ? (10^3, "μs") :
               timecenter < 10^10 ? (10^6, "ms") :
                                   (10^9, " s")
           @sprintf("%4.1f%s [%4.1f%s, %4.1f%s]",
               timecenter/tscale, tunits, timelower/tscale, tunits, timeupper/tscale, tunits)
       end

       "Five-point-stencil finite-difference approximation of the bivariate Laplacian."
       function laplacian(N)
           thediag = speye(N, N)
           offdiag = sparse(2:N, 1:(N-1), 1, N, N)
           tridiag = 2*thediag - offdiag - offdiag'
           kron(tridiag, thediag) + kron(thediag, tridiag)
       end

       A = laplacian(10^3); # matrix order is N^2 = 10^6!

       Benchmarks.@benchmarkable(nativeamd, nothing, ApproxMinimumDegree.amd(A), nothing)
       nativeamdres = Benchmarks.execute(nativeamd, 30, 30)
       Benchmarks.@benchmarkable(wrappedamd, nothing, AMD.amd(A), nothing)
       wrappedamdres = Benchmarks.execute(wrappedamd, 30, 30)
       
       println("\nApproxMinimumDegree : $(prettytimes(nativeamdres))")
       println("AMD                 : $(prettytimes(wrappedamdres))")

       end;

ApproxMinimumDegree : 633.5ms [624.5ms, 642.4ms]
AMD                 : 613.8ms [582.9ms, 644.8ms]
```
Room for optimization on this demo problem remains.

ApproxMinimumDegree.jl does not yet perform aggressive absorption nor dense row delay. As
such, ApproxMinimumDegree.jl should perform poorly relative to AMD and HSL-MC47 where
those optimizations significantly impact runtime. Neither aggressive absorption nor
dense row delay should significantly impact the preceding demo.

## Future
(Beyond a better interface and additional optimizations.)

ApproxMinimumDegree.jl was designed with two overarching features in mind. The first feature
is abstraction and separation of the common components of local preordering methods --- the
quotient graph model, the node weight priority queue, and the node weight updating scheme.

On the one hand, this abstraction and decoupling cost some peformance relative to local
preordering codes that admit tight coupling, though the second feature ---
cache-friendliness, insofar as possible --- mitigates this performance cost.

On the other hand, this abstraction and decoupling is ApproxMinimumDegree.jl's future:
Recyling local preordering methods' major infrastructure (the quotient graph and node weight
priority queue) enables rapid instantiation of high-performance local preordering methods by
leaving only the node weight updating scheme for a new local preordering method to be written.

Hence the long-term vision for ApproxMinimumDegree.jl is evolution into LocalPreorderings.jl,
a high-performance framework for experimentation with local preordering methods and a
collection of established methods.

## References

ApproxMinimumDegree.jl was written with extensive reference to the following materials
(though not associated code):

- I.S. Duff, A.M. Erisman, and J.K. Reid. "Direct methods for sparse matrices." Oxford: Clarendon press, 1986.

- I.S. Duff and J.K. Reid. "The multifrontal solution of indefinite sparse symmetric linear." ACM Transactions on Mathematical Software (TOMS) 9.3 (1983): 302-325.

- I.S. Duff and J.K. Reid. "MA27 --- a set of Fortran subroutines for solving sparse symmetric sets of linear equations." UKAEA Atomic Energy Research Establishment, 1982.

- P.R. Amestoy, T.A. Davis, and I.S. Duff. "An approximate minimum degree ordering algorithm." SIAM Journal on Matrix Analysis and Applications 17.4 (1996): 886-905.

- P.R. Amestoy, T.A. Davis, and I.S. Duff. "Algorithm 837: AMD, an approximate minimum degree ordering algorithm." ACM Transactions on Mathematical Software (TOMS) 30.3 (2004): 381-388.

- G. Reißig. "Local fill reduction techniques for sparse symmetric linear systems." Electrical Engineering 89.8 (2007): 639-652.

- C. Ashcraft. "Compressed graphs and the minimum degree algorithm." SIAM Journal on Scientific Computing 16.6 (1995): 1404-1411.

- A. George and J.W.H. Liu. "The evolution of the minimum degree ordering algorithm." Siam review 31.1 (1989): 1-19.

- A. George and J.W.H. Liu. "Computer solution of large sparse positive definite systems." Prentice-Hall Inc., Englewood Cliffs, N.J., 1981.

- A. George and J.W.H. Liu. "A fast implementation of the minimum degree algorithm using quotient graphs." ACM Transactions on Mathematical Software (TOMS) 6.3 (1980): 337-358.

- A. George and J.W.H. Liu. "A minimal storage implementation of the minimum degree algorithm." SIAM Journal on Numerical Analysis 17.2 (1980): 282-299.

- A. George and D.R. McIntyre. "On the application of the minimum degree algorithm to finite element systems." Springer Berlin Heidelberg, 1977.