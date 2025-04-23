### spline3

Simple implementation of cubic spline interpolation. This project is mostly to get an idea of the workflow utilizing unit tests (with googletest).

Cubic splines are piecewise cubic polynomials, defined with respect to subintervals. They use a local representation of the form $p(h) = a + b h + c h^2 + d h^3$, where $h$ is the local variable (relative to left sub-interval end-point). In this implementation, the spline interpolation is formulated in terms of a system for the $b$-coefficients on all intervals (i.e. the first derivatives there).

Warning: It seems to be more conventional to solve in terms of the $c$-coefficients. Please take this pitfall into account, when comparing the system with known expressions. 

Implemented boundary conditions are:
* clamped
* not-a-knot

Further notes:
- The code should work on arbitrary grids.
- No effort at optimization has been made.
