#' Different version of HRK Procedure
#'
#' See [est_arma_hrk].
#' Stage III of the Hannan-Rissanen-Kavalieris procedure is implemented here as well.
#' This function returns the best model (since the iterations might not always improve the log-likelihood value) and
#' allows for returning results at different stages of the HRK procedure.
#' \cr
#' One notable differences is that the data needs to be demeaned because
#' we use Yule-Walker estimation in the first stage to ensure stability and because
#' the implementation of stage III would otherwise be even more cumbersome.
#'
#' @inheritParams est_arma_hrk
#' @param tol_stage2,tol_stage3 Default set to `1e-3`.
#'     Maximal absolute distance of the deep parameters of adjacent iterations
#' @param maxit_stage2,maxit_stage3 Integers. Default for both stages is 5.
#' @param info Boolean. Indicates whether the slot `extra_info` of the returned list
#'     should contain a tibble with additional info about the model at different iterations,
#'     e.g., min and max absolute value of zeros, poles, value of objective function etc.
#' @param tol Small tolerance, used to check whether the mean of the supplied data matrix is indeed zero.
#'
#' @return See [est_arma_hrk()].
#'   The list contains the additional slots
#'     \describe{
#'       \item{stage_opt}{Since the best stable and miniphase model is returned, we also indicate whether this happened at stage II or stage III.}
#'       \item{info_tibble}{Tibble containing all relevant info about the outcome at different stages and iterations.}
#'     }
#'  and does not contain the slots `y.mean` (since mean is required to be zero) and `converged`.
#' @export
#'
#' @examples
#' data = BQdata_xts
#' set.seed(123)
#' for (pp in 0:2){
#'   for (qq in 0:2){
#'     if (pp + qq == 0){next}
#'     cat(paste0("p = ", pp, " q = ", qq, "\n"))
#'     tmpl = tmpl_arma_pq(m = 2, n = 2,
#'                         p = pp, q = qq)
#'     ff = est_arma_hrk3(data, tmpl = tmpl)
#'   }
#' }
#' ff$info_tibble %>% print(n=100)
est_arma_hrk3= function(y, tmpl,
                        maxit_stage2 = 5, tol_stage2 = 1e-3,
                        maxit_stage3 = 5, tol_stage3 = 1e-3,
                        info = TRUE,
                        trace = FALSE,
                        p.max = NULL, ic = c("AIC", "BIC", "max"),
                        mean_estimate = c("zero", "sample.mean", "intercept"),
                        tol = sqrt(.Machine$double.eps)) {

  # Input checks and retrievals ####

  ic = match.arg(ic)
  mean_estimate = match.arg(mean_estimate)
  stopifnot("est_arma_hrk3(): Input argument *mean_estimate* needs to be equal to *zero*." = mean_estimate == "zero",
            "est_arma_hrk3(): Mean not zero (within tolerance). Data need to be demeaned." = all(colMeans(y) < tol))

  # ***This can be deleted after checking that it isn't used anywhere
  intercept = 0

  # save original template
  tmpl_orig = tmpl
  ok = is.template(tmpl_orig, 'armamod')
  stopifnot("est_arma_hrk3(): the template has no free/deep parameters!" = tmpl_orig$n.par > 0)

  # Check dimensions
  order = unname(tmpl_orig$order)
  dim_out = order[1] # output dimension
  dim_in = order[1] # input dimension
  deg_ar = order[3] # degree of a(z)
  deg_ma = order[4] # degree of b(z)

  stopifnot("est_arma_hrk3(): Only square ARMA models are allowed." = dim_out == dim_in,
            "est_arma_hrk3(): The number of outputs must be greater than zero." = dim_out > 0,
            "est_arma_hrk3(): The model must not be constant (AR order or MA order must be non-zero)" = max(deg_ar, deg_ma) > 0)


  # Adjust template: Not all templates are allowed ####

  # The following must be satisfied
  # * a[0] = b[0]
  # * diagonal elements of a[0] (and b[0]) are equal to 1
  # * all other fixed elements are zero!
  # * sigma_L is lower triangular
  #
  # A new template is created which takes the above into account and subsequently compared to the originally supplied template

  # Linear indices for all elements in (a_0, ... , a_p, b_0, ... , b_q), written in wide matrix form and column-wise
  ix_lin = matrix(1:((dim_out^2)*(deg_ar+deg_ma+2)), nrow = dim_out, ncol = dim_out*(deg_ar+deg_ma+2))

  # New template for AR and MA parameters
  ab_tmpl = matrix(0, nrow = dim_out, ncol = dim_out*(deg_ar+deg_ma+2))

  stopifnot("est_arma_hrk3(): The original template and the newly created one are not of the same length." = length(ab_tmpl) == (length(tmpl_orig$h) - dim_out^2))

  # Find the linear indices of the free elements (the respective row in the affine trafomatrix H must not be zero if the element is free)
  free = apply((tmpl_orig$H[1:((dim_out^2)*(deg_ar+deg_ma+2)), , drop = FALSE] != 0), MARGIN = 1, any)
  ab_tmpl[free] = NA_real_

  # Set diagonal elements of a[0] equal to one
  diag(ab_tmpl) = 1

  # Fix b[0] here as intermediate step. Below it will be made equal to a[0])
  ab_tmpl[, (dim_out*(deg_ar+1)+1):(dim_out*(deg_ar+2))] = 0

  # Create new template (including matrices (h,H) describing the restrictions)
  tmpl = model2template(armamod(structure(ab_tmpl, order = c(dim_out, dim_out, deg_ar, deg_ma), class = c('lmfd', 'ratm')),
                                sigma_L = matrix(NA_real_, nrow = dim_out, ncol = dim_out)),
                        sigma_L = 'chol')

  # Take the constraint b[0] = a[0] into account: Copy the (h,H) part pertaining to a[0] to the part pertaining to b[0]
  ix_lin_a = ix_lin[, 1:dim_out]
  ix_lin_b = ix_lin[, (dim_out*(deg_ar+1)+1):(dim_out*(deg_ar+2))]
  tmpl$h[ix_lin_b]   = tmpl$h[ix_lin_a]
  tmpl$H[ix_lin_b, ] = tmpl$H[ix_lin_a, ]

  # Compare the original and the constructed template
  if (!isTRUE(all.equal(cbind(tmpl_orig$h, tmpl_orig$H), cbind(tmpl$h, tmpl$H)))) {
    warning('The restrictions of the supplied template has been adapted to satisfy
            diag(a_0) = 1, a_0 = b_0,
            fixed elements need to be fixed to zero,
            sigma_L needs to be lower-triangular.')
  }

  # Data ####

  # Check data input
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('Input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }

  # Sample size
  stopifnot("est_arma_hrk3(): Data input argument *y* needs to have non-trivial length" = length(y) > 0,
            "est_arma_hrk3(): Data input argument *y* is not compatible with the template *tmpl*." = ncol(y) == dim_out,
            "est_arma_hrk3(): Data input argument *y* contains NAs/NaNs/Infs." = all(is.finite(y)))
  n_obs = nrow(y)


  # Stage I: Long AR Regression ####

  # estimate the noise by a estimating long AR model
  AR_estimate = TRUE

  # Set maximal AR order for selecting the optimal one
  if (is.null(p.max)) {

    # Sometimes, p.max is given as the max of 12 and the term below.
    # For 10*log10(n_obs) to be smaller than 12, the sample size *n_obs* would need to be 15 or smaller
    p.max = min(10*log10(n_obs), floor((n_obs-1)/(dim_out+1)))
  }
  p.max = as.integer(p.max)[1]

  stopifnot("est_arma_hrk3(): Stage I: Maximal AR order *p.max* is negative. Can only happen if sample size is smaller than 10." = p.max >= 0,
            "est_arma_hrk3(): Stage I: Sample size *n_obs* is not sufficient for the desired maximum AR order *p.max*" = (n_obs-p.max) >= (dim_out*p.max + intercept))

  # Set the model selection criterion
  if (ic == 'max') penalty = -1
  if (ic == 'AIC') penalty = 2/n_obs
  if (ic == 'BIC') penalty = log(n_obs)/n_obs

  # Estimate the AR model with Yule-Walker (for stability reasons)
  out = est_ar(y, p.max = p.max, penalty = penalty, method = "y")

  # Recalculate residuals and obtain innovation covariance
  e = solve_inverse_de(sys = out$model$sys, y)$u
  sigma1 = crossprod(e)/nrow(e)

  # Poles of a(z)
  z_a = try(out$model$sys$a %>% zeroes(print_message = FALSE) %>% roots_as_list()) ### WS ###
  z_a_abs = try(z_a %>% purrr::map_dbl(abs))

  # ___*Extrainfo: Stage 1: Residuals, cov and its trace, zeros of a(z) and their absolute value  ####
  if(info){
    info_tibble = tibble::tibble(stage = 1, iteration = 1, flipped = FALSE,
                                 min_pole = min(z_a_abs), max_pole = max(z_a_abs),
                                 min_zero = Inf, max_zero = Inf,
                                 trace = sum(diag(sigma1)), log_lik = Inf, lndetSigma = out$stats[out$p,3])
  } else {
    # Provide a default value
    info_tibble = tibble::tibble()
  }

  # Stage II: Construct inputs ####

  # __Observations and their Lags as (block) Toeplitz matrix ####
  # The matrix Y below won't change any more!
  R = array(NA_real_, dim=c(1, dim_out, deg_ar+1))
  R[ , , 1] = y[1, ]
  C = t(y)
  dim(C) = c(1, dim_out, n_obs)
  Y = btoeplitz(R, C)

  # __Construct regressor matrix (-y_t, -y_{t-1}, ... , -y_{t-p}, e_t, ... , e_{t-q})
  # Note that y_t and e_t are also included here to allow for contemporaneous dependence modeling.
  # Identification/exclusion restrictions will of course be taken care of!
  # The matrix E will change during the iterations!
  # The dimensions of R stay the way they are here. Errors are written into the R array in the iterative stages
  R = array(NA_real_, dim=c(1, dim_out, deg_ma+1))
  E = matrix(NA_real_, nrow = n_obs, ncol = (deg_ma+1)*dim_out)
  YE = cbind(-Y, E)

  # Print additional info
  if (trace) {
    cat('HRK estimation of ARMA model: ',
        'm =', dim_out, ', n_obs=', n_obs, ', p=', deg_ar, ', q=', deg_ma, '\n', sep = '')
    if (AR_estimate) {
      cat('initial AR estimate of noise p.max=', p.max,
          ', p=', out$p, ', ll=', (-1/2)*(dim_out*(log(2*pi) + 1) + out$stats[out$p+1, 'lndetSigma']),
          '\n', sep = '')
    }
    cat('iter |th - th0|  n.val      MSE       ll \n')
  }

  # Stage II: Iterating ####

  # Initialize parameters for stage II
  iter_stage2 = 1
  iterate_stage2 = TRUE

  # Optimization outcomes
  th = numeric(tmpl$n.par)
  ll_opt = Inf
  ll0 = dim_out * (log(2*pi) + 1)
  model_opt = NULL

  while (iterate_stage2 && iter_stage2 <= maxit_stage2) {

    # Initial value for deep parameters (for comparison)
    th0 = th

    # Create regression matrix of contemporaneous and lagged residuals
    # Note that R is of dim=c(1, dim_out, deg_ma+1), and that "R" overwrites the contents of "C" in btoeplitz(),
    # therefore only the diagonal block is relevant, the rest in "R" is NAs
    R[, , 1] = e[1, ]
    C = t(e)
    dim(C) = c(1, dim_out, n_obs)
    E = btoeplitz(R, C)

    # __Combine all data: LHS and RHS ####

    # Write errors "E" from just above into the regression matrix YE (which also contains contemporaneous values)
    YE[, ((deg_ar+1)*dim_out+1):((deg_ar+deg_ma+2)*dim_out)] = E
    # Compare HD page 296, Stage IIa, point (ii): In the reversed echelon form,
    # the last variables in the vector of outputs have to be regressed
    # on contemporaneous values of (\hat{\varepsilon_t} - y_t), such that the coefficients of a[0] and b[0] coincide, plus the other variables.
    # Here, it is more general in the sense that we are not constrained to the reversed Kronecker echelon form
    # but only to a certain set of affine restrictions on the system and noise parameters
    #
    # To state the same in formulas, we have
    # a[0] y[t] + a[1] y[t-1] + ... = a[0] e[t] + b[1] e[t-1] + ...
    # y[t] = a*[0] (e[t] - y[t]) + a[1] (-y[t-1]) + ... + b[1] e[t-1] + ... + e[t]
    # where a*[0] = (a[0] - I), which is zero if a[0] is the identity matrix
    YE[ , 1:dim_out] = e - y

    # Obtain valid observations and initialize regression residuals for univariate regressions
    ix_no_na = apply(!is.na(cbind(y, YE)), MARGIN = 1, FUN = all)
    n.valid = sum(ix_no_na)

    # Initialize new residuals
    u = y
    u[!ix_no_na, ] = NA_real_

    # Note that ab_tmpl has zeros for b[0]!
    ab = ab_tmpl

    # __Row-wise regressions ####
    for (ix_var in (1:dim_out)) {

      ix_free = which(is.na(ab_tmpl[ix_var, ]))
      if (length(ix_free) > 0) {

        stopifnot("est_arma_hrk3(): Stage II regression: Valid sample size is not sufficient for the desired ARMA model." = sum(ix_no_na) >= (length(ix_free)+intercept))

        # __HERE IS THE REGRESSION ####
        # YE contains dim_out*(deg_ar + deg_ma + 2) columns (also contemporaneous outputs and inputs)
        # The fact that b[0] = a[0] and diag(a[0]) = 1 is taken care of
        out = stats::lsfit(YE[ix_no_na, ix_free], y[ix_no_na, ix_var], intercept = intercept)
        ab[ix_var, ix_free] = out$coefficients
        u[ix_no_na, ix_var] = out$residuals
      }
    }

    # set b[0] = a[0]
    # Remember that the first block column in the regression matrix YE is (e-y)
    ab[, (dim_out*(deg_ar+1)+1):(dim_out*(deg_ar+2))] = ab[, 1:dim_out]

    # estimated ARMA system
    sys = structure(ab,
                    order = c(dim_out, dim_out, deg_ar, deg_ma),
                    class = c('lmfd', 'ratm'))

    # Covariance of residuals of ARX regression
    sigma2 = crossprod(u[ix_no_na, , drop = FALSE])/n.valid
    sigma2_L = t(chol(sigma2))

    # Zeros of b(z) and a(z)
    z_b = sys$b %>% zeroes(print_message = FALSE) %>% roots_as_list()  ### WS ###
    z_b_abs = z_b %>% purrr::map_dbl(abs)
    z_a = sys$a %>% zeroes(print_message = FALSE) %>% roots_as_list()  ### WS ###
    z_a_abs = z_a %>% purrr::map_dbl(abs)

    # ___*Extrainfo: Stage 2 zeroes of b(z) and a(z) ####
    if (info){
      info_tibble =
        dplyr::bind_rows(info_tibble,
                         tibble::tibble(stage = 2, iteration = iter_stage2, flipped = FALSE,
                                        min_pole = min(z_a_abs), max_pole = max(z_a_abs),
                                        min_zero = min(z_b_abs), max_zero = max(z_b_abs),
                                        trace = sum(diag(sigma2)), log_lik = -1/2*(ll0 + log(det(sigma2))), lndetSigma = log(det(sigma2))))
    }

    # __Flipping zeros of b(z) outside ####
    b = sys$b

    if (any(z_b_abs < 1)){
      # Associated error process should have identity matrix as covariance
      b = b %r% sigma2_L

      # Flip all zeros inside the unit circle outside
      b_flipped = reflect_zeroes(b, z_b[which(z_b_abs < 1)] %>% unlist())

      # Normalize b(z) such that b(0) = I: Postmultiply poly b(z) with b(0)^{-1}
      b_flipped_at_zero = zvalue(b_flipped, z = 0)
      b2 = b_flipped %r% solve(b_flipped_at_zero) %r% ab[, 1:dim_out, drop = FALSE]

      # Reassign
      sys = lmfd(a = sys$a, b = b2)
    }

    # Obtain inputs as b^{-1}(z) a(z) y_t: This explodes if b(z) is non-miniphase!
    e = solve_inverse_de(sys, y)$u

    # __noise covariance and log-lik value ####
    sigma2f = crossprod(e[ix_no_na, , drop=FALSE])/n.valid
    sigma2f_L = try(t(chol(sigma2f)))

    if (inherits(sigma2f_L, 'try-error')) {
      stop('Estimated noise covariance matrix is not positive definite.\n',
           'This is probably due to a non minimum phase estimate of the ARMA system.')
    } else {
      det_sigma2f_L = prod(diag(sigma2f_L))
      ll = ll0 + 2*log(det_sigma2f_L)
      error = ''
    }

    # estimated model
    model = armamod(sys, sigma_L = sigma2f_L)

    # estimated parameters
    th = extract_theta(model, tmpl, on_error = "stop")

    # Print additional info
    if (trace) {
      cat(sprintf('%4.0f %10.3f %6.0f %8.3f %8.3f',
                  iter_stage2, max(abs(th-th0), na.rm = TRUE), n.valid, sum(diag(sigma2f)),
                  (-1/2)*ll), error, '\n')
    }

    # Zeros of b(z) and a(z)
    z_b = sys$b %>% zeroes(print_message = FALSE) %>% roots_as_list() ### WS ####
    z_b_abs = z_b %>% purrr::map_dbl(abs)

    # ___*Extrainfo: Stage 2 flipped: Residuals, cov and its trace ####
    if (info){
      info_tibble =
        dplyr::bind_rows(info_tibble,
                         tibble::tibble(stage = 2, iteration = iter_stage2, flipped = TRUE,
                                        min_pole = min(z_a_abs), max_pole = max(z_a_abs),
                                        min_zero = min(z_b_abs), max_zero = max(z_b_abs),
                                        trace = sum(diag(sigma2)), log_lik = -1/2*(ll0+2*log(det_sigma2f_L)), lndetSigma = 2*log(det_sigma2f_L)))
    }

    # __ Optimal model: Is it stable and better than the previous iteration? ####
    if (all(z_a_abs > 1) &&  ll < ll_opt){
      model_opt = model
      th_opt = th
      residuals_opt = e
      sigma_opt = sigma2f
      ll_opt = ll
      stage_opt = 2
      iter_opt = iter_stage2
    }

    iter_stage2 = iter_stage2 + 1
    iterate_stage2 = iterate_stage2 && (iter_stage2 <= maxit_stage2) && (max(abs(th - th0)) > tol_stage2)
  }

  if ((trace) && (max(abs(th - th0)) <= tol_stage2)) {
    cat('algorithm converged in stage 2\n')
  }

  # Stage III ####

  iter_stage3 = 1
  iterate_stage3 = TRUE
  if(trace){ cat("\nStage 3: \n") }

  # In case maxit_stage3 == 0, we need default values for sigma3f
  sigma3f = sigma2f

  while (iterate_stage3 && iter_stage3 <= maxit_stage3){

    th0 = th

    # _Define LHS ####

    # Solve b(z) (\eta_t, \xi_t) = (y_t, e_t)
    eta_lhs = solve_de(sys = lmfd(a = sys$b, b = diag(dim_out)), y)$y
    xi_lhs = solve_de(sys = lmfd(a = sys$b, b = diag(dim_out)), e)$y

    # Obtain result
    lhs = e + eta_lhs - xi_lhs

    # Vectorize
    lhs = c(t(lhs))

    # _Define RHS ####

    # Allocate
    rhs_eta = matrix(NA_real_, nrow = n_obs*dim_out, ncol = dim_out^2)
    rhs_xi = matrix(NA_real_, nrow = n_obs*dim_out, ncol = dim_out^2)


    # The RHS contains b(z)^-1 [ (y_t', ... , y_{t-p}') %x% I_n ] stacked one above the other (same for stage II residuals)
    for(ix_pos in 1:dim_out){
      for (ix_var in 1:dim_out) {

        # Input for solve_de()
        input = matrix(0, nrow = n_obs, ncol = dim_out)

        # Component ix_var of the observations is at position ix_pos
        input[, ix_pos] = y[, ix_var]

        new_regr_eta = solve_de(sys = lmfd(a = sys$b, b = diag(dim_out)),
                                u = input)$y
        rhs_eta[, (ix_var - 1)*dim_out + ix_pos] = -c(t(new_regr_eta))

        # Same same for residuals
        input[, ix_pos] = e[, ix_var]

        new_regr_xi = solve_de(sys = lmfd(a = sys$b, b = diag(dim_out)),
                               u = input)$y
        rhs_xi[, (ix_var-1)*dim_out + ix_pos] = c(t(new_regr_xi))
      }
    }

    # Lags for rhs_eta and stacking in the right way
    R = array(NA_real_, dim=c(dim_out, dim_out^2, deg_ar+1))
    R[ , , 1] = rhs_eta[1:dim_out, , drop = FALSE]
    C = t(rhs_eta)
    dim(C) = c(dim_out^2, dim_out, n_obs)
    C = aperm(C, perm = c(2,1,3))
    RHS_ETA = btoeplitz(R, C)

    # Lags for rhs_xi and stacking in the right way
    R = array(NA_real_, dim=c(dim_out, dim_out^2, deg_ma+1))
    R[ , , 1] = rhs_xi[1:dim_out, , drop = FALSE]
    C = t(rhs_xi)
    dim(C) = c(dim_out^2, dim_out, n_obs)
    C = aperm(C, perm = c(2,1,3))
    RHS_XI = btoeplitz(R, C)

    # Binding the columns to complete the regressor matrix on the right hand side
    rhs = cbind(RHS_ETA, RHS_XI)

    # # Obtain lagged values
    # # iseq() because deg_ar could be zero
    # for (ix_lag_ar in iseq(1, deg_ar)){
    # rhs[, (dim_out^2)*ix_lag_ar + (1:dim_out^2)] =
    #     xts::lag.xts(rhs[, (1:(dim_out^2)), drop = FALSE], k = ix_lag_ar*dim_out)
    # }
    #
    # # deg_ma is greater than one
    # for (ix_lag_ma in 1:(deg_ma)) {
    #   rhs[, (dim_out^2*(deg_ar + 1)) + (dim_out^2)*ix_lag_ma + (1:dim_out^2)] =
    #     xts::lag.xts(rhs[,  (dim_out^2*(deg_ar + 1)) + (1:dim_out^2), drop = FALSE], k = ix_lag_ma*dim_out)
    # }

    # Contemporaneous regressors: Take a[0] = b[0] into account
    rhs[, 1:(dim_out^2)] = rhs[, 1:(dim_out^2)] + rhs[, dim_out^2*(deg_ar+1) + 1:(dim_out^2)]

    # __Regression ####
    ab = ab_tmpl
    ab_vec = c(ab)
    ix_free = which(is.na(ab_vec))
    ix_no_na_vec = apply(!is.na(cbind(lhs, rhs)[, ix_free, drop = FALSE]), MARGIN = 1, FUN = all)
    n_valid = sum(ix_no_na_vec)
    u = lhs
    u[!ix_no_na_vec] = NA_real_

    stopifnot("Some parameters must be free" = length(ix_free) > 0,
              "Number of valid observations must be larger than number of regressors" = n_valid > length(ix_free) + intercept,
              "mean_estimate = 'intercept' not yet implemented for Stage 3. Please use 'sample.mean' instead." = intercept == FALSE)

    lsfit_out_stage3 = stats::lsfit(rhs[ix_no_na_vec, ix_free, drop = FALSE],
                                    lhs[ix_no_na_vec],
                                    intercept = intercept)

    # Update ARMA system (lmfd)
    ab[ix_free] = lsfit_out_stage3$coefficients

    # Remember that the first block column in the regression matrix YE is (e-y)
    ab[, (dim_out*(deg_ar+1)+1):(dim_out*(deg_ar+2))] = ab[, 1:dim_out]

    # estimated ARMA system
    sys = structure(ab,
                    order = c(dim_out, dim_out, deg_ar, deg_ma),
                    class = c('lmfd', 'ratm'))

    # New residuals for i-th variable (was initialized with the outputs)
    u[ix_no_na_vec] = lsfit_out_stage3$residuals
    u = t(matrix(u, nrow = dim_out, ncol = n_obs))
    ix_no_na = apply(!is.na(u), MARGIN = 1, FUN = all)

    sigma3 = crossprod(u[ix_no_na, , drop = FALSE])/sum(ix_no_na)
    sigma3_L = t(chol(sigma3))

    z_b = sys$b %>% zeroes(print_message = FALSE) %>% roots_as_list() ### WS ###
    z_b_abs = z_b %>% purrr::map_dbl(abs)

    z_a = sys$a %>% zeroes(print_message = FALSE) %>% roots_as_list() ### WS ###
    z_a_abs = z_a %>% purrr::map_dbl(abs)

    # ___*Extrainfo: Stage 3 zeroes of b(z) ####
    if(info){
      info_tibble = dplyr::bind_rows(info_tibble,
                                     tibble::tibble(stage = 3, iteration = iter_stage3, flipped = FALSE,
                                                    min_pole = min(z_a_abs), max_pole = max(z_a_abs),
                                                    min_zero = min(z_b_abs), max_zero = max(z_b_abs),
                                                    trace = sum(diag(sigma3)), log_lik = -1/2*(ll0+2*log(prod(diag(sigma3_L)))), lndetSigma = 2*log(prod(diag(sigma3_L))) ))

    }

    # __Flip unstable zeros of b(z) outside if necessary ####
    b = sys$b

    if (any(z_b_abs < 1)){
      # Associated error process should have identity matrix as covariance
      b = b %r% sigma3_L

      # Flip all zeros inside the unit circle outside
      b_flipped = reflect_zeroes(b, z_b[which(z_b_abs < 1)] %>% unlist())

      # Normalize b(z) such that b(0) = I: Postmultiply poly b(z) with b(0)^{-1} and a_0
      b_flipped_at_zero = zvalue(b_flipped, z = 0)
      b3 = b_flipped %r% solve(b_flipped_at_zero) %r% ab[, 1:dim_out, drop = FALSE]

      # Reassign
      sys = lmfd(a = sys$a, b = b3)
    }

    # __Calculate residuals with solve_de() ####
    # (No intercept here)
    e = solve_inverse_de(sys, y)$u

    # __Calculate cov-mat and likelihood value ####
    sigma3f = crossprod(e[ix_no_na, , drop=FALSE])/sum(ix_no_na)
    sigma3f_L = try(t(chol(sigma3f)))
    if (inherits(sigma3f_L, 'try-error')) {
      stop('Estimated noise covariance matrix is not positive definite.\n',
           'This is probably due to a non minimum phase estimate of the ARMA system.')
    } else {
      det_sigma3f_L = prod(diag(sigma3f_L))
      ll = ll0 + 2*log(det_sigma3f_L)
      error = ''
    }

    # Zeros of Stage 3: b(z) and a(z)
    z_b = sys$b %>% zeroes(print_message = FALSE) %>% roots_as_list() ### WS ###
    z_b_abs = z_b %>% purrr::map_dbl(abs)

    # ___*Extrainfo: Stage 3: Flipped residuals, their cov-mat, and trace ####
    if (info){
      info_tibble =
        dplyr::bind_rows(info_tibble,
                         tibble::tibble(stage = 3, iteration = iter_stage3, flipped = TRUE,
                                        min_pole = min(z_a_abs), max_pole = max(z_a_abs),
                                        min_zero = min(z_b_abs), max_zero = max(z_b_abs),
                                        trace = sum(diag(sigma3f)),
                                        log_lik = -1/2*(ll0+2*log(prod(diag(sigma3f_L)))),
                                        lndetSigma = 2*log(prod(diag(sigma3f_L)))))
    }

    # __Update model, extract deep parameters ####
    model = armamod(sys, sigma_L = sigma3f_L)
    th = extract_theta(model, tmpl, on_error = "stop")

    # Print additional info
    if (trace) {
      cat(sprintf('%4.0f %10.3f %6.0f %8.3f %8.3f',
                  iter_stage3, max(abs(th-th0), na.rm = TRUE), n.valid, sum(diag(sigma3f)),
                  (-1/2)*ll), error, '\n')
    }

    # __ Optimal model: Is it stable and better than the previous iteration? ####
    if (all(z_a_abs > 1) &&  ll < ll_opt){
      model_opt = model
      th_opt = th
      residuals_opt = e
      sigma_opt = sigma3f
      ll_opt = ll
      stage_opt = 3
      iter_opt = iter_stage3
    }

    iter_stage3 = iter_stage3 + 1
    iterate_stage3 = iterate_stage3 && (iter_stage3 <= maxit_stage3) && (max(abs(th - th0)) > tol_stage3)
  }

  stopifnot("estimate_arma_hrk(): No iteration in stage II or III resulted in a stable model!" = !is.null(model_opt))

  return(list(model = model_opt, th = th_opt, tmpl = tmpl,
              residuals = residuals_opt, sigma = sigma_opt, n.valid = n.valid, ll = (-1/2)*ll_opt,
              stage = stage_opt, iter = iter_opt,
              info_tibble = info_tibble) )
}




#' Hannan, Rissanen, Kavalieris estimation procedure
#'
#' `r lifecycle::badge("deprecated")`
#' Estimate (V)ARMA models with the Hannan, Rissanen Kavalieris procedure, see
#' e.g. \insertCite{Hannan.Rissanen82}{RLDM} and \insertCite{HannanKavalierisMackisack1986}{RLDM}.
#'
#' The main idea of the HRK procedure is as follows. If we have given estimates,
#' \eqn{e_t}{e[t]} say, of the disturbances,
#' then the ARMA model is estimated from the equation
#' \deqn{y_t = -a^*_0(y_t+e_t) - a_1 y_{t-1} - \cdots - a_p y_{t-p} +
#'                 b_1 e_{t-1} + \cdots + b_q e_{t-q} + v_{t-1}}{
#'       y[t] = -a^*[0] (y[t]+e[t]) - a[1] y[t-1] - \dots - a[p] y[t-p] +
#'                 b[1] e[t-1] + \dots + b[q] e[t-q] + v[t]}
#' where \eqn{a^*_0}{a^*[0]} is obtained from \eqn{a_0 = b_0}{a[0]=b[0]} by setting all
#' diagonal elements equal to zero.
#' The entries in the parameter matrices \eqn{a_i}{a[i]} and \eqn{b_i}{b[i]} are
#' either treated as fixed (and equal to zero)
#' or "free". Now the above regression is estimated "componentwise", i.e. for
#' each component of \eqn{y_t}{y[t]} the corresponding
#' "free" parameters are estimated by OLS.
#'
#' Given the parameter estimates one computes new estimates for the disturbances,
#' by recursively solving the ARMA system,
#' see [solve_inverse_de()]. The sample variance of these residuals is
#' used as an estimate of the
#' noise covariance matrix \eqn{\Sigma}.
#'
#' This procedure may be iterated: use the "new" estimates for the disturbances to (re) estimate
#' the ARMA parameters and to (re) estimate the disturbances, ...
#'
#' The parameters `maxit` and `tol` control this iterative scheme. The
#' iterations are stopped after at most `maxit` iterations or when there is
#' only a "small" change of the estimates. To be more precise, if `th`,
#' `th0` denote the vector of parameter estimates in the actual round and
#' the previous round, then the procedure stops if `max(abs(th-th0)) <=
#' tol`.
#'
#' Note that in general there is no guarantee that this iterative scheme
#' converges or that the estimates are improved by iterating.
#'
#' The user may supply his "own" (initial) estimates `e` of the disturbances.
#' If the parameter `e` is missing (or `NULL`) then the procedure `est_arma_hrk`
#' computes estimates of the disturbances by fitting a "long" AR model to the data. To this end
#' the procedure simply calls [est_ar_ols()] with the respective paramaters
#' `p.max` (which controls the maximum possible AR order), `ic` (which controls the
#' information criterion used to select the order of the AR model) and `mean_estimate`
#' (which tells `est_ar_ols` how to estimate the mean \eqn{\mu}).
#' The default for the maximum order `p.max` is
#' \deqn{\max(12, 10\log_{10}(N), (N-1)/(m+1))}{max(12, 10 log[10](N), (N-1)/(m+1))}
#'
#' The procedure supports three options for the estimation of the mean
#' \eqn{\mu = \mathbf{E} y_t}{\mu = E y[t]}. For `mean_estimate="zero"` the
#' procedure sets the (estimate of the) mean equal to zero. For
#' `mean_estimate="sample.mean"` the procedure simply uses the sample mean
#' of `y` as an estimate. Third option `mean_estimate="intercept"`
#' uses an intercept in the above regression(s) and computes the estimate of the
#' mean correspondingly. Note that this fails if the estimated AR polynomial has
#' a unit root, i.e. if \deqn{\det \hat{a}(1) = 0.}{det a(1) = 0.}
#'
#' There is no guarantee that the HRK algorithm returns a stable and minimum
#' phase ARMA model. In particular, if the estimated model is *not* minimum
#' phase then the recursive computation of the residuals often yields useless
#' results and correspondingly the cholesky decomposition of the sample variance
#' of the residuals (which is used as estimate of the noise covariance matrix
#' \eqn{\Sigma}) fails. In this case the procedure stops with an error message.
#'
#' @param y sample, i.e. an \eqn{(N,m)} dimensional matrix,
#'   or a "time series" object (i.e. `as.matrix(y)` should return an
#'   \eqn{(N,m)}-dimensional numeric matrix). Missing values (`NA`, `NaN` and
#'   `Inf`) are **not** supported.
#' @param e (initial) estimate of the disturbances \eqn{u_t}{u[t]}.
#'   If non `NULL` then `e` has to be an \eqn{(N,m)} dimensional matrix,
#'   or a "time series" object (i.e an object which may be coerced to an \eqn{(N,m)}
#'   dimensional matrix with `as.matrix(e)`).
#'   \cr
#'   If `NULL` then the procedure computes an estimate of the disturbances by fitting
#'   a "long" AR model to the data, see [est_ar_ols()].
#'   \cr
#'   The matrix `e` may contain missing values (`NA`, `NaN` and
#'   `Inf`). Note that e.g. `est_ar_ols` returns residuals where the first
#'   \eqn{p} (here \eqn{p} refers to the order of the fitted AR model) values are missing.
#' @param tmpl a model template, see [model structures()]. Note that only the case
#'   is implemented, where
#'   \eqn{a_0=b_0}{a[0]=b[0]} holds, the diagonal entries of
#'   \eqn{a_0=b_0}{a[0]=b[0]} are equal to one and all other fixed elements
#'   are equal to zero. Furthermore the square root `sigma_L` of the noise covariance matrix
#'   is asssumed to be a lower triangular matrix without any further restrictions.
#'   \cr
#'   The given template is coerced to a template of this kind. If the given template does not
#'   comply to these restrictions, then a warning message is issued.
#' @param maxit (integer) maximum number of iterations
#' @param tol (numeric) tolerance level
#' @param trace (boolean) if `trace=TRUE`, then some tracing information on the
#'              iterations is printed.
#' @param p.max (integer or `NULL`) Maximum order of the candidate AR models.
#'     For the default choice see below.
#' @param ic (character string) Which information criterion shall be used to find the
#'     optimal AR order. Note that `ic="max"` means that an AR(p) model
#'     with `p=p.max` is fitted. Default is `ic="AIC"`.
#' @param mean_estimate Character string giving the method used to estimate the mean \eqn{\mu}.
#'     Default is `mean_estimate = "sample.mean"`.
#'     See the details below.
#'
#' @return List with components
#' \item{model}{the estimated (V)ARMA model (i.e. an [armamod()] object).}
#' \item{th}{vector with the (free) parameters of the estimated (V)ARMA model.}
#' \item{tmpl}{the (coerced) model template.}
#' \item{y.mean}{estimate of the mean \eqn{\mu}.}
#' \item{residuals}{the residuals of the model, computed with [solve_inverse_de()].}
#' \item{sigma}{the sample variance \eqn{S} of the residuals, i.e. an estimate of the
#'              noise covariance matrix \eqn{\Sigma}.}
#' \item{n.valid}{number of "valid" observations, i.e. observations where all needed
#'                lagged values \eqn{y_{t-i}}{y[t-i]} and \eqn{e_{t-i}}{e[t-i]}
#'                are availiable. For an ARMA(p,q) model this implies
#'                that the number of valid observations is less than or equal
#'                to `n.obs -max(p,q)`.}
#' \item{ll}{Gaussian log likelihood:
#'           \deqn{(-1/2)(m \ln(2\pi) + m + \ln\det(S))}{(-N/2)(m ln(2\pi) + m + ln det(S))}
#'           where \eqn{S} denotes the sample variance of the residuals.}
#' \item{iter}{number of iterations.}
#' \item{converged}{(boolean) indicates whether the algorithm converged.}
#'
#' @export
#'
#' @references
#' \insertRef{Hannan.Rissanen82}{RLDM}
#'
#' \insertRef{HannanKavalierisMackisack1986}{RLDM}
#'
#' @examples
#' # in order to get reproducible results
#' set.seed(4321)
#'
#' # generate a random VARMA(p=2,q=1) model with m=2 outputs #####################
#' tmpl = tmpl_arma_pq(m = 2, n = 2, p = 2, q = 1)
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#'
#' # generate a sample with 200 observations
#' data = sim(model, n.obs = 200, n.burn_in = 100)
#'
#' # estimate model with HRK
#' # note: we are cheating here and use the true disturbances!
#' out = est_arma_hrk(data$y, data$u, tmpl)
#' print(out$model)
#' # ll() returns the same logLik value. However, we have to demean the data
#' all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
#'                      'conditional', skip = 2))
#'
#' # estimate the model with HRK
#' # use the residuals of a long AR model as estimates for the noise
#' out = est_arma_hrk(data$y, e = NULL, tmpl,
#'                    trace = TRUE, maxit = 10, mean_estimate = 'zero')
#' print(out$model)
#' # ll() returns the same logLik value. However, we have to demean the data
#' all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
#'                      'conditional', skip = 2))
#'
#' # Generate a random Model in echelon form model (m = 3) #######################
#' tmpl = tmpl_arma_echelon(nu = c(1,1,1))
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#'
#' # generate a sample with 200 observations
#' data = sim(model, n.obs = 200, n.burn_in = 100)
#' # add mean value(s)
#' data$y = data$y + matrix(1:3, nrow = 200, ncol = 3, byrow = TRUE)
#'
#' # estimate model with HRK
#' # note: we are cheating here and use the true disturbances!
#' out = est_arma_hrk(data$y, data$u, tmpl,
#'                    trace = FALSE, maxit = 1, mean_estimate = 'sample.mean')
#' print(out$y.mean)
#' print(out$model)
#' # ll() returns the same logLik value. However, we have to demean the data
#' all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
#'                      'conditional', skip = 1))
#'
#' # estimate the model with HRK
#' # use the residuals of a long AR model as estimates for the noise
#' out = est_arma_hrk(data$y, e = NULL, tmpl,
#'                    maxit = 10, mean_estimate = 'intercept')
#' print(out$y.mean)
#' print(out$model)
#' # ll() returns the same logLik value. However, we have to demean the data
#' all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
#'                      'conditional', skip = 1))
#'
#' # We may also use this procedure to estimate AR models #####################
#' # where some coefficients are fixed = 0
#' a = dbind(d = 3, diag(2), array(NA_real_, dim = c(2,2,2)))
#' a[1,2,] = 0 # all coefficient matrices are lower triangular, i.e.
#' # y[2t] does not Granger cause y[1t]
#' tmpl = model2template(armamod(sys = lmfd(a=a),
#'                               sigma_L = matrix(NA_real_, nrow = 2, ncol = 2)),
#'                       sigma_L = 'chol')
#' model = r_model(template = tmpl, bpoles = 1, bzeroes = 1, sd = 0.25)
#' diag(model$sigma_L) = 1 # scale the diagonal entries of sigma_L
#' print(model)
#'
#' # generate a sample with 200 observations
#' data = sim(model, n.obs = 200, n.burn_in = 100)
#'
#' # estimate model with HRK
#' out = est_arma_hrk(data$y, NULL, tmpl,
#'                    trace = FALSE, maxit = 1, mean_estimate = 'zero')
#' print(out$y.mean)
#' print(out$model)
#' # ll() returns the same logLik value. However, we have to demean the data
#' all.equal(out$ll, ll(out$model, scale(data$y, center = out$y.mean, scale = FALSE),
#'                      'conditional', skip = 2))
#'
#' # reset the "seed"
#' set.seed(NULL)
est_arma_hrk = function(y, e = NULL, tmpl,
                        maxit = 1, tol = 1e-3, trace = TRUE,
                        p.max = NULL, ic = c("AIC", "BIC", "max"),
                        mean_estimate = c("sample.mean", "intercept", "zero")) {

  ic = match.arg(ic)
  mean_estimate = match.arg(mean_estimate)
  intercept = (mean_estimate == 'intercept')

  # check template "tmpl"
  tmpl0 = tmpl # save original template
  ok = is.template(tmpl0, 'armamod')

  order = unname(tmpl0$order)
  if ( (order[1] != order[2]) || (order[1] < 1)  || (max(order[c(3,4)]) < 1) ) {
    stop('only square, ARMA models with m > 0 and (p+q) > 0 are implemented!')
  }
  if (tmpl0$n.par == 0) stop('the template has no free/deep parameters!')

  m = order[1] # output dimension
  p = order[3] # degree of a(z)
  q = order[4] # degree of b(z)

  # Template ####
  # Coerce the template to a simplified template according to
  #  diagonal elements of a[0] and b[0] are equal to 1
  #  all other fixed elements are zero!
  #  a[0] = b[0]
  #  sigma_L is lower triangular
  i = matrix(1:((m^2)*(p+q+2)), nrow = m, ncol = m*(p+q+2))
  ab_tmpl = matrix(0, nrow = m, ncol = m*(p+q+2)) # AR,MA parameters

  if (length(ab_tmpl) != (length(tmpl0$h) - m^2)) stop('this should not happen!')
  free = apply((tmpl0$H[1:((m^2)*(p+q+2)), , drop = FALSE] != 0), MARGIN = 1, any)
  ab_tmpl[free] = NA_real_
  # diagonal elements of a[0] are equal to one
  diag(ab_tmpl) = 1
  # b[0] = a[0], therefore we fix b[0] here
  ab_tmpl[, (m*(p+1)+1):(m*(p+2))] = 0
  # print(ab_tmpl)

  # create a new template, according to this structure:
  model = armamod(structure(ab_tmpl, order = c(m,m,p,q), class = c('lmfd', 'ratm')),
                  sigma_L = matrix(NA_real_, nrow = m, ncol = m))
  tmpl = model2template(model, sigma_L = 'chol')
  # take the constraint b[0] = a[0] into account
  ia = i[,1:m]
  ib = i[,(m*(p+1)+1):(m*(p+2))]
  tmpl$h[ib]   = tmpl$h[ia]
  tmpl$H[ib, ] = tmpl$H[ia, ]

  # compare the given and the constructed template
  if (!isTRUE(all.equal(cbind(tmpl0$h, tmpl0$H), cbind(tmpl$h, tmpl$H)))) {
    warning('the given template has been adapted to the restrictions!')
  }

  # Data ####
  # coerce data object(s) to matrices
  y = try(as.matrix(y))
  if ( inherits(y, 'try-error') || (!is.numeric(y)) || (!is.matrix(y)) ) {
    stop('input "y" must be a data object which may be coerced to a matrix with "as.matrix(y)')
  }
  if (length(y) == 0) stop('"y" contains no data')
  n.obs = nrow(y)  # sample size
  if (ncol(y) != m) stop('"y" is not compatible with the template "tmpl"!')
  if (any(!is.finite(y))) stop('"y" contains NAs/NaNs/Infs!')


  # __Residual estimates ####
  AR_estimate = FALSE

  # Estimate long AR regression if no residuals are supplied
  if (is.null(e)) {
    if ( (q == 0) && (!any(is.na(ab_tmpl[, 1:m]))) && all(ab_tmpl[ ,1:m]== diag(m)) ) {
      # pure AR model, we don't need initial estimate of the noise
      e = matrix(0, nrow = n.obs, ncol = m)
    } else {
      # estimate the noise by a estimating long AR model
      AR_estimate = TRUE

      if (is.null(p.max)) {
        p.max = max(1, min(10*log10(n.obs)/m, floor((n.obs-1)/(m+1))))
      }
      p.max = as.integer(p.max)[1]
      if (p.max < 0) stop('p.max must be a non negative integer!')
      if ( (n.obs-p.max) < (m*p.max + intercept) ) {
        stop('sample size N is not sufficient for the desired maximum AR order "p.max"')
      }

      if (ic == 'max') penalty = -1
      if (ic == 'AIC') penalty = 2/n.obs
      if (ic == 'BIC') penalty = log(n.obs)/n.obs

      out = est_ar_ols(y, p.max = p.max, penalty = penalty, mean_estimate = mean_estimate)
      e = out$residuals
    }
  }

  # Check residuals (input or obtained from long AR regression)
  e = try(as.matrix(e))
  if ( inherits(e, 'try-error') || (!is.numeric(e)) || (!is.matrix(e)) ) {
    stop('input "e" must be a data object which may be coerced to a matrix with "as.matrix(e)')
  }
  if ((ncol(e) != m) || (nrow(e) != n.obs)) stop('matrix of residuals "e" is not compatible!')

  # __Mean adjust observations ####
  y.mean = double(m)
  if (mean_estimate == 'sample.mean') {
    y.mean = colMeans(y, na.rm = TRUE)  # sample mean
    y = y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE)
  }

  # __Lagged data: Toeplitz matrices ####
  # construct (block) Toeplitz matrices with the lagged data
  R = array(NA_real_, dim=c(1, m, p+1))
  R[ , , 1] = y[1, ]
  C = t(y)
  dim(C) = c(1, m, n.obs)
  Y = btoeplitz(R, C)
  # print(Y[1:10,])

  # __Concatenate data matrices in the right format
  R = array(NA_real_, dim=c(1, m, q+1)) # Errors are written into this array in the iterative stages
  E = matrix(NA_real_, nrow = n.obs, ncol = (q+1)*m)
  YE = cbind(-Y, E)

  # Print additional info
  if (trace) {
    cat('HRK estimation of ARMA model: ',
        'm=', m, ', n.obs=', n.obs, ', p=', p, ', q=', q, '\n', sep = '')
    if (AR_estimate) {
      cat('initial AR estimate of noise p.max=', p.max,
          ', p=', out$p, ', ll=', (-1/2)*(m*(log(2*pi) + 1) + out$stats[out$p+1, 'lndetSigma']),
          '\n', sep = '')
    }
    cat('iter |th - th0|  n.val      MSE       ll \n')
  }

  # Iterative Stages: ####

  # Initialize counters etc
  iter = 1
  iterate = TRUE
  th = numeric(tmpl$n.par)
  ll0 = m * (log(2*pi) + 1)

  #
  while (iterate) {
    # Initial value for deep parameters
    th0 = th

    # Create regression matrix of lagged residuals
    R[, , 1] = e[1, ] # note that R is of dim=c(1, m, q+1)
    C = t(e)
    dim(C) = c(1, m, n.obs)
    E = btoeplitz(R, C)
        # print(E[1:10,])

    # combine all data
    # a[0] y[t] + a[1] y[t-1] + ... = a[0] e[t] + b[1] e[t-1] + ...
    # y[t] = a*[0] (e[t] - y[t]) + a[1] (-y[t-1]) + ... + b[1] e[t-1] + ... + e[t]
    # Note that a*[0] is zero if a[0] is the identity matrix
    YE[, ((p+1)*m+1):((p+q+2)*m)] = E # write errors from just above into the big regression matrix (with observations)
    YE[ , 1:m] = e - y

    # Obtain valid observations and intialise regression residuals for univariate regressions
    v = apply(!is.na(cbind(y, YE)), MARGIN = 1, FUN=all)
    n.valid = sum(v)
    u = y
    u[!v, ] = NA_real_

    ab = ab_tmpl
    # __Row-wise regressions ####
    for (i in (1:m)) {
      k = which(is.na(ab_tmpl[i, ]))
      # cat('i', i, 'dim(YE)', dim(YE), 'length(v)', length(v), 'k', k,'\n')
      # print(ab_tmpl[i, ])
      if (length(k) > 0) {
        if (sum(v) < (length(k)+intercept)) {
          stop('sample size "N" is not sufficient for the desired ARMA model')
        }
        out = stats::lsfit(YE[v, k], y[v, i], intercept = intercept) # note that YE contains m(p+q+2) columns!!! The fact that b[0] = a[0] and diag(a[0]) = 1 is taken care of!
        if (intercept) {
          y.mean[i] = out$coefficients[1]
          ab[i, k] = out$coefficients[-1]
        } else {
          ab[i, k] = out$coefficients
        }
        u[v, i] = out$residuals
      }
    }
    # set b[0] = a[0]
    ab[, (m*(p+1)+1):(m*(p+2))] = ab[, 1:m]
    # mean estimate
    if (intercept) {
      a1 = ab[, 1:(m*(p+1)), drop = FALSE]
      dim(a1) = c(m,m,p+1)
      a1 = apply(a1, MARGIN = c(1,2), FUN = sum)
      # y.mean0 = y.mean
      y.mean = solve(a1, y.mean)
      # print(cbind(y.mean0, a1, solve(a1), y.mean))
    }
    # estimated ARMA system
    sys = structure(ab, order = c(m, m, p, q), class = c('lmfd', 'ratm'))
      # (re)estimate noise, given this estimate of the ARMA system
      if (intercept) {
        e = solve_inverse_de(sys, y - matrix(y.mean, nrow = n.obs, ncol = m, byrow = TRUE))$u
      } else {
        e = solve_inverse_de(sys, y)$u
      }

    # noise covariance
    sigma = ( t(e[v, , drop=FALSE]) %*% e[v, , drop=FALSE] )/n.valid
    # print(sigma)
    # print(eigen(sigma)$values)
    sigma_L = try(t(chol(sigma)))
    if (inherits(sigma_L, 'try-error')) {
      stop('Estimated noise covariance matrix is not p.d.\n',
           'This is probably due to a non minimum phase estimate of the ARMA system.')
      # ll = NA_real_
      # sigma_L = diag(x = sqrt(diag(sigma)), nrow = m, ncol = m)
      # iterate = FALSE
      # error = 'sigma is not p.d.'
    } else {
      det_sigma_L = prod(diag(sigma_L))
      ll = ll0 + 2*log(det_sigma_L)
      error = ''
    }

    # estimated model
    model = armamod(sys, sigma_L = sigma_L)

    # estimated parameters
    th = extract_theta(model, tmpl)

    # Print additional info
    if (trace) {
      cat(sprintf('%4.0f %10.3f %6.0f %8.3f %8.3f',
                  iter, max(abs(th-th0), na.rm = TRUE), n.valid, sum(diag(sigma)),
                  (-1/2)*ll), error, '\n')
      # cat(max(abs(th - th_TRUE)), '\n')
      # cat(y.mean - y.mean_TRUE, '\n')
    }

    iter = iter + 1
    iterate = iterate && (iter <= maxit) && (max(abs(th - th0)) > tol)
  }

  if ((trace) && (max(abs(th - th0)) <= tol)) {
    cat('algorithm converged\n')
  }

  return( list(model = model, th = th, tmpl = tmpl, y.mean = y.mean,
               residuals = e, sigma = sigma, n.valid = n.valid, ll = (-1/2)*ll,
               iter = iter - 1, converged = (max(abs(th - th0)) <= tol) ) )
}


