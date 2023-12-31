% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{outputs_ARMA_cpp}
\alias{outputs_ARMA_cpp}
\title{Outputs of an ARMA systems}
\usage{
outputs_ARMA_cpp(A1, B, t0, u, y)
}
\arguments{
\item{A1}{\eqn{(m, mp)} matrix \eqn{-a_0^{-1}(a_p,...,a_1)}{-a[0]^{-1}(a[p],...,a[1])}.}

\item{B}{\eqn{(m, n(q+1))} matrix \eqn{a_0^{-1}(b_0,...,b_q}{a[0]^{-1}(b[0],...,b[q])}.}

\item{t0}{integer, start iteration at t = t0.}

\item{u}{\eqn{(n, N)} matrix with the inputs \eqn{(u_1,...,u_N}{(u[1],...,u[N])}.}

\item{y}{\eqn{(m, N)} matrix with the outputs \eqn{(y_1,...,y_N}{(y[1],...,y[N])}.}
}
\value{
This RcppArmadillo routine returns \code{NULL} but \bold{overwrites}
the input argument \code{y} with the computed outputs!
}
\description{
This internal helper function computes the outputs of an ARMA system
\deqn{a_0 y_t + a_1 y_{t-1} + \cdots + a_p y_{t-p} = b_0 u_t + \cdots + b_q u_{t-q}}{
      a[0] y[t] + a[1] y[t-1] + ... + a[p] y[t-p] = b[0] u[t] + ... + b[q] u[t-q]}
}
\details{
Values \eqn{y_t}{y[t]}, \eqn{u_t}{u[t]} for \eqn{t\leq 0}{t\le 0} are implicitly set to be zero.
However, by starting the iteration with some \eqn{t_0>1}{t0>1} we can enforce non-zero
initial values.
}
\note{
Use this procedure with care!
\itemize{
\item The procedure does \bold{not} check the input arguments. We require \eqn{m > 0},
\eqn{p \geq 0}{p \ge 0}, \eqn{n(q+1) \geq 0}{n(q+1)\ge 0} and
\eqn{1 \leq t_0 \leq N}{1 \le t_0  \le N}.
\item The procedure \bold{overwrites} the input argument \code{y}.
\item The data matrices are organized columnwise (to avoid memory shuffling)!
\item Note also the non standard representation of the coefficient matrices.
}
}
\examples{
# generate a random ARMA(2,1) model (3 outputs, 2 inputs)
p = 2
q = 1
m = 3
n = 2
model = test_armamod(dim = c(m, n), degrees = c(p,q), digits = 2)
A = unclass(model$sys$a)
a0 = A[,,1]
A1 = -A[,,(p+1):2]
dim(A1) = c(m, m*p)
A1 = solve(a0, A1)
B = unclass(model$sys$b)
dim(B) = c(m, n*(q+1))
B = solve(a0, B)

# generate random noise sequence (sample size N = 10)
n.obs = 10
u = matrix(rnorm(n.obs*n), nrow = n, ncol = n.obs)
print(u)

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# call outputs_ARMA_cpp()
outputs_ARMA_cpp(A1, B, t0 = 2, u, y) # start with t>=2
print(u)
print(y)  # y is overwritten with the computed outputs
}
\seealso{
\code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
\code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
\code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
}
