% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{residuals_STSP_cpp}
\alias{residuals_STSP_cpp}
\title{Residuals of a statespace system}
\usage{
residuals_STSP_cpp(A, B, C, D, y, a, u, dPI, dU)
}
\arguments{
\item{A}{\eqn{(s,s)} matrix.}

\item{B}{\eqn{(s,m)} matrix.}

\item{C}{\eqn{(m,s)} matrix.}

\item{D}{\eqn{(m,m)} matrix, must be regular.}

\item{y}{\eqn{(m,N)} matrix with the outputs: \eqn{(y_1,y_2,\ldots,y_N)}{(y[1],y[2],...,y[N])}.}

\item{a}{\eqn{(s,N+1)} matrix. This matrix is overwritten with the (computed) states:
\eqn{(a_1,a_2,\ldots,a_N,a_{N+1})}{(a[1],a[2],\ldots,a[N],a[N+1])}.
On input \code{a[,1]} must hold the initial state \eqn{a_1}{a[1]}.}

\item{u}{\eqn{(m,N)} matrix. This matrix is overwritten with (computed) residuals:
\eqn{(u_1,u_2,\ldots,u_N)}{(u[1],u[2],...,u[N])}.}

\item{dPI}{\eqn{((m+s)^2,K)} matrix.}

\item{dU}{\eqn{(mN,K)} matrix or \eqn{(0,0)} matrix. This matrix is overwritten with the
directional derivatives of the residuals. However, if
the matrix is empty then no derivatives are computed.}
}
\value{
This RcppArmadillo implementation returns \code{NULL} but \bold{overwrites} the input
arguments \code{a}, \code{u} and \code{dU}!
}
\description{
This internal helper function computes the residuals
(and the directional derivatives of the residuals)
for a statespace system of the form
\deqn{a_{t+1} = A a_t + B u_t, \; y_t = C a_t + D u_t}{
      a[t+1] = A a[t] + B u[t] and y[t] = C a[t] + D u[t]}
The system must be square and non-empty, i.e. \eqn{m=n>0}.
}
\note{
Use this procedure with care!
\itemize{
\item The procedure does \bold{not} check the input arguments.
\item The procedure \bold{overwrites} the input arguments
\code{a}, \code{u} and \code{dU}.
\item The data matrices are organized columnwise (to avoid memory shuffling)!
}
}
\examples{
# generate a random statespace model (3 outputs, 3 inputs and 4 states)
m = 2
s = 3
model = test_stspmod(dim = c(m, m), s = s, digits = 2)

# generate random noise sequence (sample size N = 10)
n.obs = 10
u = matrix(rnorm(n.obs*m), nrow = m, ncol = n.obs)

# generate matrix for the state sequence
a = matrix(0, nrow = s, ncol = n.obs+1)
a[,1] = rnorm(s) # random initial state a[1]

# generate matrix for the outputs
y = matrix(0, nrow = m, ncol = n.obs)

# compute outputs and states
outputs_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D, u, a, y)

# recompute the states and disturbances/residuals from the given outputs:
uu = u + 0 # "deep copy" of the disturbances
aa = a + 0 # and the states
aa[, 2:(n.obs+1)] = 0 # clear all states a[t], t > 1
residuals_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D,
                   y, aa, uu, diag(0), diag(0))
all.equal(u, uu) # check
all.equal(a, aa) # check

# compute directional derivatives of residuals
dPI = diag((m+s)^2)
dU = matrix(0, nrow = n.obs*m, ncol = ncol(dPI))
residuals_STSP_cpp(model$sys$A, model$sys$B, model$sys$C, model$sys$D,
                   y, aa, uu, dPI, dU)

# check the directional derivatives
eps = 1e-8
dU_num = matrix(0, nrow = m*n.obs, ncol = (m+s)^2)
dPI = matrix(0, nrow = (m+s), ncol = (m+s))
for (k in (1:((m+s)^2))) {
  dPI[] = 0
  dPI[k] = eps
  uu = u + 0 # "deep copy" of the disturbances
  aa = a + 0 # and the states
  residuals_STSP_cpp(model$sys$A + dPI[1:s,1:s],
                     model$sys$B + dPI[1:s,(s+1):(s+m)],
                     model$sys$C + dPI[(s+1):(s+m),1:s],
                     model$sys$D + dPI[(s+1):(s+m),(s+1):(s+m)],
                     y, aa, uu, diag(0), diag(0))
  dU_num[, k] = c(uu - u )/eps  # num. approx. of the derivative in direction "dPI"
}
# relative error of the numerical approximation
junk = (abs(dU)+abs(dU_num))
junk[junk == 0] = 1
2*abs(dU_num - dU)/junk
}
\seealso{
\code{\link{outputs_ARMA_cpp}}, \code{\link{residuals_ARMA_cpp}},  \code{\link{cll_theta_ARMA_cpp}},
\code{\link{outputs_STSP_cpp}}, \code{\link{residuals_STSP_cpp}},  \code{\link{cll_theta_STSP_cpp}} and
\code{\link{solve_de}}, \code{\link{solve_inverse_de}} and \code{\link{ll}}.
}
