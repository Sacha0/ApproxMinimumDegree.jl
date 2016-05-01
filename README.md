# ApproxMinimumDegree.jl

ApproxMinimumDegree.jl is a native-Julia implementation of the approximate minimum external
degree algorithm. This is a prototype; prefer [AMD](http://faculty.cse.tamu.edu/davis/suitesparse.html)
([AMD.jl](https://github.com/JuliaOptimizers/AMD.jl)) or [HSL-MC47](http://www.hsl.rl.ac.uk/catalogue/mc47.html)
for anything beyond experimentation.

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