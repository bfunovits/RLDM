% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aa_rldm_classes.R
\name{armamod}
\alias{armamod}
\title{Constructor for LMFD (ARMA) Models}
\usage{
armamod(sys, sigma_L = NULL, names = NULL, label = NULL)
}
\arguments{
\item{sys}{\code{\link[rationalmatrices:lmfd]{rationalmatrices::lmfd()}} or \code{\link[rationalmatrices:rmfd]{rationalmatrices::rmfd()}} object}

\item{sigma_L}{Left-factor of noise covariance,
i.e. the covariance \eqn{\sigma} is obtained as \code{sigma_L * t(sigma_L)}.
If \code{sigma_L} is a vector of dimension \eqn{n}, where \eqn{n} is the input dimension, only the diagonal elements are parametrized.
If it is a vector of dimension \eqn{n^2}, then the elements of \code{sigma_L} are filled column by column.}

\item{names}{optional vector of character strings}

\item{label}{optional character string}
}
\value{
Object of class \code{armamod}.
}
\description{
A left-matrix fraction description (LMFD) plus parameterisation of noise covariance.
}
\details{
In Hannan, Deistler (2012, page 7), RMFDs are also called dynamic adjustment forms.
Internally, MFDs are lists with slots \code{sys}, \code{sigma_L}, \code{names}, \code{label}.
}
\examples{
x = armamod(sys = lmfd(c(1, 0.5), 1), sigma_L = diag(1))
x
}
