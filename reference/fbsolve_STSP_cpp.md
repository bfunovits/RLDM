# Forward-backward solution of statespace systems

**DEPRECATED**? This internal helper function computes the outputs of
an, in general **unstable**, statespace system \$\$a\_{t+1} = A a_t + B
u_t, \\ y_t = C a_t + D u_t\$\$ by forward and backward recursion. The
procedure assumes that the state transition matrix \\A\\ is block upper
triangular, where the upper block \\A\_{11}\\ is stable (i.e. all
eigenvalues have moduli less than one) and the lower block \\A\_{22}\\
is unstable (i.e. all eigenvalues have moduli larger than one). This
function is mainly used in the routine
[`innovation_form`](https://bfunovits.github.io/RLDM/reference/innovation_form.md).

## Usage

``` r
fbsolve_STSP_cpp(A, B, C, D, u, au, as, y)
```

## Arguments

- A:

  \\(s, s)\\ matrix.

- B:

  \\(s, n)\\ matrix.

- C:

  \\(m, s)\\ matrix.

- D:

  \\(m, n)\\ matrix.

- u:

  \\(n, N)\\ matrix with the inputs \\(u_1,...,u_N\\.

- au:

  \\(su,N+1)\\ matrix. This matrix is **overwritten** with the
  (computed) states of the unstable part of the system.
  \\(a\_{u1},a\_{u2},\ldots,a\_{uN},a\_{u,N+1})\\. On input `au[,N+1]`
  must hold the "initial" state \\a\_{u,N+1}\\.

- as:

  \\(ss,N+1)\\ matrix. This matrix is **overwritten** with the
  (computed) states of the stable part of the system.
  \\(a\_{s1},a\_{s2},\ldots,a\_{sN},a\_{s,N+1})\\. On input `as[,1]`
  must hold the "initial" state \\a\_{s1}\\.

- y:

  \\(m,N)\\ matrix. This matrix is **overwritten** with (computed)
  outputs: \\(y_1,y_2,\ldots,y_N)\\.

## Value

This RcppArmadillo routine returns `NULL` but **overwrites** the input
argument `y`, `au` and `as` with the computed outputs and states!

## Note

Use this procedure with care!

- The procedure does **not** check the input arguments. We require \\m
  \> 0\\, \\n \> 0\\. Furthermore it is assumed that the state
  transition matrix \\A\\ is block upper triangular, as explained above.

- The procedure **overwrites** the input arguments `y`, `as` and `au`.

- The data matrices are organized columnwise (to avoid memory
  shuffling)!

## See also

[`outputs_STSP_cpp`](https://bfunovits.github.io/RLDM/reference/outputs_STSP_cpp.md)
and
[`innovation_form`](https://bfunovits.github.io/RLDM/reference/innovation_form.md).
