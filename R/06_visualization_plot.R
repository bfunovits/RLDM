# plot.____ methods ####################################################################
#

#' Plot Methods
#'
#' Plot methods for impulse response functions ([impresp()] objects),
#' autocovariance functions ([autocov()] objects),
#' frequency response functions ([freqresp()] objects) and
#' spectral densities ([spectrald()] objects).
#'
#' If `x` is an \eqn{(m,n)} dimensional object then the plot is divided into an
#' \eqn{(m,n)} array of subfigures. In each each of the subfigures the respective
#' \eqn{(i,j)}-th element of the object `x` is displayed. The methods allow
#' simultaneous plotting of several objects, by passing a list `x_list`
#' with additional objects to the procedure. In the following we assume that
#' `x_list` contains \eqn{k-1} objects, i.e. in total \eqn{k} objects
#' are to be plotted.
#'
#' The parameter "`xlim`" determines the x-limits of the subfigures:
#' \cr
#' `xlim='global'` uses the same x-limits for all subfigures (the limits are
#' determined from the data). In the case of frequency response or spectral density
#' objects one may also pass a two-dimensional vector `xlim = c(x1,x2)` to
#' the plot method. In this case all subfigures use these values as common x-limits.
#' \cr
#' `xlim='column'` means that all sub figures in a "column" have
#' the same x-limits. The limits are determined from data.
#' Finally `xlim='subfigure'` means that each subfigure gets its own
#' x-limits (determined from data).
#'
#' Quite analogously the parameter `ylim` determines the limits for the y-axes.
#' (Just replace `'column'` by `'row'`).
#'
#' The plot methods have quite a number of optional design parameters. In most cases
#' these parameters are interpreted as follows.
#' `NA` values mean that the methods use some suitable defaults.
#' E.g. the labels for the x- and y-axis are chosen according to the class of the
#' object(s) to be plotted and the parameter "`which`".
#' `NULL` values for the (optional) parameters mean that the respective graphic element is omitted.
#' E.g. `subfigure_main = NULL` skips the titles for the subfigures.
#'
#' The titles for the (m,n) subfigures are determined by the parameter `subfigure_main`.
#' One may pass an \eqn{(m,n)} character matrix or a scalar character (string) to the procedure.
#' If `subfigure_main` is a scalar (character string) then the procedures creates
#' the respective titles by replacing the "place holders" `i_` and `j_` with the
#' respective row and column number. See the examples below.
#'
#' The "style" parameters `col, type, ..., bg.points` determine the appearance of the
#' "lines" for each of the `k` objects. (If necessary these values are "recycled".)
#' See also [graphics::lines()] and [graphics::points()] for a
#' more detailed explanation of these parameters.
#'
#' If more than one object is to be plotted (the optional parameter `x_list` is not empty)
#' then a suitable legend may be added with the parameters `legend` and
#' `legend_args`. Note that `legend` should be a character (or expression) vector of
#' length `k`.
#'
#' The parameter "`which`" determines what to plot in the case of
#' "`freqresp`" or "`spectrald`" objects.
#' \describe{
#' \item{gain,modulus}{plot the moduli `abs(x[i,j])` versus frequencies.}
#' \item{phase}{plot the arguments `Arg(x[i,j])` versus frequencies.}
#' \item{nyquist}{plot the imaginary part `Im(x[i,j])` versus the real part `Re(x[i,j])`.}
#' \item{coherence}{plot the *coherence*. This plot is somewhat special.
#'       The "coherence matrix" is symmetric and and the diagonal entries are
#'       are equal to one (for all frequencies). Therefore the entries on and below the
#'       diagonal contain no additional information.
#'       For this reason the subfigures above the diagonal display the coherence,
#'       the subfigures below show the "scaled arguments" `Arg(x[i,j])/(2*pi)+0.5`
#'       and the subfigures on the diagonal display the scaled auto spectra of the
#'       \eqn{m} component processes.}
#' }
#'
#' These plot methods use the internal helper function [rationalmatrices::plot_3D()].
#'
#' @param x  [impresp()], [autocov()], [freqresp()] or
#'           [spectrald()]  object.
#' @param x_list (optional) list of additional objects (of the same class as "`x`").
#' @param xlim,ylim determine the axis limits of the subfigures.
#'          E.g. `xlim = 'column'` means that all subfigures in a column
#'          use the same x-axis limits. Analogously `y = 'row'` implies
#'          that the subfigures in a row share the same limits for the y-axis.
#'          \cr
#'          For the "`freqresp`" and "`spectrald`" plots the parameter
#'          `xlim` may also be a numeric 2-dimensional vector `xlim = c(x1,x2)`.
#'          In this case all sub-figures use the given limits for the x-axis.
#'          Furthermore the limits for the y-axis are computed based on the
#'          corresponding data subset. This option may be used to "zoom"
#'          into a certain range of frequencies.
#'
#' @param sampling_rate (number) sampling rate.
#' @param unit (character string) time or frequency unit.
#' @param which (character string) what to plot. This parameter is only used
#'              for plotting frequency response objects and spectral densities.
#'              See details below.
#' @param log a character string which contains "x" if the x axis is to be logarithmic,
#'            "y" if the y axis is to be logarithmic and "xy" or "yx" if both axes
#'            are to be logarithmic.  This parameter is only used
#'            for plotting frequency response objects and spectral densities.
#'            Note that a logarithmic y-axis only makes sense when
#'            plotting the moduli of the frequency response or the spectral density.
#' @inheritParams rationalmatrices::plot_3D
#' @param ... not used.
#'
#' @return The plot methods return (invisibly) a "[closure()]", `subfig` say, which may be used to
#'         add additional graphic elements to the subfigures. The call `opar = subfig(i,j)`
#'         creates a new (sub) plot at the (i,j)-th position with suitable margins and
#'         axis limits. See the examples below.
#' @importFrom graphics par plot rect text
#' @export
#'
#' @rdname plot
#' @name plot methods
#'
#' @examples
#' set.seed(1995)  # set seed to get reproducible results
#' n.obs = 2^8
#' m = 2
#' s = 3
#'
#' # generate a random, stable and minimum phase state space model
#' # for a bivariate process (x[t], y[t])
#' model = test_stspmod(dim = c(m,m), s = s, bpoles = 1, bzeroes = 1)
#' model$names = c('x[t]', 'y[t]')
#'
#' # simulate data
#' data = sim(model, n.obs = n.obs)
#'
#' #### plot impulse response
#' # overlay three different "orthogonalization" schemes
#' plot(impresp(model), list(impresp(model, H = 'eigen'),
#'                           impresp(model, H = 'chol')),
#'      legend = c('none','chol','eigen'),
#'      legend_args = list(title = 'orthogonalization', fill = NA, border = NA, bty = 'n'),
#'      style = 'colored', ylim = 'subfig', xlab = NA)
#'
#' #### plot partial autocorrelation function
#' # overlay with the corresponding sample partial ACF
#'
#' par(lend = 1) # in order to get a "barplot"
#' plot(autocov(data$y, lag.max = 12, type = 'partial'),
#'      list(autocov(model, lag.max = 12, type = 'partial')),
#'      subfigure_main = 'delta[i_*j_](k)', parse_subfigure_main = TRUE,
#'      style = 'bw', type = c('h','l'), pch = 19, lwd = c(15,2),
#'      legend = c('sample', 'true'))
#' par(lend = 0) # reset 'lend=0'
#'
#' # frequency response of the model
#' n.f = 2^11
#' frr = freqresp(model, n.f = n.f)
#'
#' #### plot "gain"
#' subfig = plot(frr, which = 'gain',
#'               sampling_rate = 60, unit = 'Hz',
#'               ylim = 'row', log = 'y',
#'               subfigure_main = 'k[i_*j_](lambda)', parse_subfigure_main = TRUE)
#'
#' # mark the frequencies with the max gain!
#' junk = unclass(frr$frr)
#' i_max = apply(Mod(junk), MARGIN = c(1,2), FUN = which.max)
#' f_max = matrix(60*((0:(n.f-1))/n.f)[i_max], nrow = 2, ncol = 2)
#' for (i in (1:2)) {
#'   for (j in (1:2)) {
#'     subfig(i,j)
#'     abline(v = f_max[i,j], col = 'steelblue')
#'   }
#' }
#'
#' #### create a "Nyquist" plot of the frequency response
#' plot(frr, which = 'nyquist',
#'      xlim = 'subfig', ylim = 'subfig',
#'      subfigure_main = 'k[i_*j_](lambda)', parse_subfigure_main = TRUE)
#'
#'
#' # compute spectral density
#' spd = spectrald(model, n.f = 256)
#'
#' #### plot the coherence
#' # the subfigure above the diagonal shows the coherenec between
#' # the two component proceses x[t] and y[t]
#' # the sub figures on the diagonal show the scaled autospectra
#' # of the two component processes x[t] and y[t].
#' # and the subfigure below the diagonal shows the
#' # phase/argument of the cross spectral density between the
#' # two component processes x[t] and y[t]
#' plot(spd, sampling_rate = 60, unit="Hz",
#'      main = expression(spectral~density~~Gamma[i*j] == kappa[i*j]*exp(i*Phi[i*j])),
#'      which = 'coherence',
#'      style = 'bw')
#'
#' # periodogram
#' per = spectrald(data$y)
#'
#' # smoothed periodogram
#' sacf = autocov(data$y, lag.max = floor(sqrt(n.obs)))
#' per2 = spectrald(sacf, n.f = 256)
#'
#'
#' #### make a plot of the absolute value of the spectral density,
#' # of the periodogram and the smoothed periogram.
#' # with a logarithmic y-axis
#' # skip zero frequency, since the periodogram is zero at lambda=0
#' plot(spd, list(per, per2), sampling_rate = 12, unit = '/year', which = 'modulus',
#'      log = 'y', xlim = c(1/n.obs, 0.5) * 12,  # skip zero frequency
#'      legend = c('true','periodogram', 'smoothed per.'),
#'      legend_args = list(bty = 'n', col = NA, lty = NA, pch = NA, lwd = 4),
#'      style = 'colored', ylim = 'subfig',
#'      subfigure_main = 'kappa[i_*j_] == group("|", Gamma[i_*j_], "|")',
#'      parse_subfigure_main = TRUE,
#'      col = c('red', 'black', 'blue'), type = 'o', lty = c(1,0,1),
#'      pch = c(NA, 19, NA), cex.points = 0.1)
#'
#' set.seed(NULL) # reset seed
plot.impresp = function(x, x_list = NULL,
                        xlim = c('global','column','subfig'), ylim = c('row','subfig','global'),
                        main = NA, xlab = NA, ylab = NULL,
                        subfigure_main = NA, parse_subfigure_main = FALSE,
                        style = c('gray', 'bw', 'bw2', 'colored'),
                        col = NA, type = 'l', lty = 'solid', lwd = 1,
                        pch = 16, cex.points = 1, bg.points = 'black',
                        legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  xlim = match.arg(xlim)
  ylim = match.arg(ylim)

  irf = unclass(x$irf)
  m0 = dim(irf)[1]
  n0 = dim(irf)[2]
  lag.max  = dim(irf)[3] - 1
  if ((m0*n0*(lag.max+1)) == 0) stop('empty impulse response')

  sigma = x$sigma_L
  sigma = sigma %*% t(sigma)
  orth = isTRUE(all.equal(sigma, diag(n0)))

  x = list(0:lag.max)
  y = list(irf)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!inherits(x_list[[i]],'impresp')) stop('argument "x_list" must be a list of "impresp" objects')
      irf = unclass(x_list[[i]]$irf)
      m = dim(irf)[1]
      n = dim(irf)[2]
      lag.max  = dim(irf)[3] - 1
      if ( (m != m0) || (n != n0)) stop('impulse response objects are not compatible')
      if ((lag.max+1) == 0) stop('empty impulse response')

      y = c(y, list(irf))
      x = c(x, list(0:lag.max))
    }
  }
  k = length(x)

  if ( (!is.null(main)) && any(is.na(main)) ) {
    if (orth) {
      main = 'orthogonalized impulse response'
    } else {
      main = 'impulse response'
    }
  }
  if ( (!is.null(xlab)) && any(is.na(xlab)) ) xlab = 'lag'
  if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = NULL

  if ( (!is.null(subfigure_main)) && any(is.na(subfigure_main)) ) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }

  subfigure = plot_3D(x, y,
                      xlim = xlim, ylim = ylim,
                      main = main, xlab = xlab, ylab = ylab,
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      legend = legend, legend_args = legend_args)
  return(invisible(subfigure))
}


#' @rdname plot
#' @export
plot.autocov = function(x, x_list = NULL,
                        xlim = c('global','column','subfig'), ylim = c('row','subfig','global'),
                        main = NA, xlab = NA, ylab = NULL,
                        subfigure_main = NA, parse_subfigure_main = FALSE,
                        style = c('gray', 'bw', 'bw2', 'colored'),
                        col = NA, type = 'l', lty = 'solid', lwd = 1,
                        pch = 16, cex.points = 1, bg.points = 'black',
                        legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  xlim = match.arg(xlim)
  ylim = match.arg(ylim)

  acf = unclass(x$acf)
  m0 = dim(acf)[1]
  lag.max  = dim(acf)[3] - 1
  if ((m0*(lag.max+1)) == 0) stop('empty impulse response')
  acf_type0 = x$type

  x = list(0:lag.max)
  y = list(acf)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!inherits(x_list[[i]],'autocov')) stop('argument "x_list" must be a list of "autocov" objects')
      acf = unclass(x_list[[i]]$acf)
      acf_type = x_list[[i]]$type
      if ( (!is.na(acf_type0)) && (acf_type != acf_type0) ) {
        acf_type0 = NA # mixed ACF types
      }
      m = dim(acf)[1]
      lag.max  = dim(acf)[3] - 1
      if ( (m != m0) ) stop('autocov objects are not compatible')
      if ((lag.max+1) == 0) stop('empty autocov object')

      x = c(x, list(0:lag.max))
      y = c(y, list(acf))
    }
  }
  k = length(x)

  if ( (!is.null(main)) && any(is.na(main)) ) {
    if (is.na(acf_type0)) {
      main = 'mixed partial autocorrelation-, autocorrelation-, autocovariance- functions'
    } else {
      main = switch(acf_type0,
                    covariance = 'autocovariance function',
                    correlation = 'autocorrelation function',
                    partial = 'partial autocorrelation function')

    }
  }
  if ( (!is.null(xlab)) && any(is.na(xlab)) ) xlab = 'lag'
  if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = NULL
  if ( (!is.null(subfigure_main)) && any(is.na(subfigure_main)) ) {
    if (m0 > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }

  subfigure = plot_3D(x, y,
                      xlim = xlim, ylim = ylim,
                      main = main, xlab = xlab, ylab = ylab,
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      cex.points = 1, bg.points = 'gray',
                      legend = legend, legend_args = legend_args)
  return(invisible(subfigure))
}


#' @rdname plot
#' @export
plot.freqresp = function(x, x_list = NULL, sampling_rate = 1, unit = '',
                         which = c('gain','phase','nyquist'),
                         xlim = NA, ylim = NA, log = '',
                         main = NA, xlab = NA, ylab = NA,
                         subfigure_main = NA, parse_subfigure_main = FALSE,
                         style = c('gray', 'bw', 'bw2', 'colored'),
                         col = NA, type = 'l', lty = 'solid', lwd = 1,
                         pch = 16, cex.points = 1, bg.points = 'black',
                         legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  which = match.arg(which)

  frr = unclass(x$frr)
  attr(frr, 'z') = NULL
  m0 = dim(frr)[1]
  n0 = dim(frr)[2]
  n.f = dim(frr)[3]
  if ((n0*m0*n.f) == 0) stop('empty frequency response')

  x = list( (0:(n.f-1))/n.f )
  y = list(frr)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!inherits(x_list[[i]],'freqresp')) stop('argument "x_list" must be a list of "freqresp" objects')

      frr = unclass(x_list[[i]]$frr)
      attr(frr, 'z') = NULL
      m = dim(frr)[1]
      n = dim(frr)[2]
      n.f = dim(frr)[3]
      if ( (m != m0) || (n != n0)) stop('freqresp objects are not compatible')
      if (n.f == 0) stop('empty frequency response')

      x = c(x, list( (0:(n.f-1))/n.f ))
      y = c(y, list(frr))
    }
  }
  k = length(x)

  if (which == 'gain') {
    for (i in (1:k)) {
      x[[i]][x[[i]] > 0.5] = x[[i]][x[[i]] > 0.5] - 1
      y[[i]] = abs(y[[i]])
      o = order(x[[i]])
      x[[i]] = sampling_rate*x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if ( (!is.null(main)) && any(is.na(main)) ) main = 'frequency response'
    if ( (!is.null(xlab)) && any(is.na(xlab)) ) {
      if (sampling_rate == 1) {
        xlab = 'frequency'
      } else {
        xlab = paste('frequency (sampling rate = ', sampling_rate, unit, ')', sep = '')
      }
    }
    if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = 'gain'
    if (any(is.na(xlim))) {
      xlim = c(0, 0.5)*sampling_rate
    }
    if (any(is.na(ylim))) {
      ylim = 'row'
    }
  }

  if (which == 'phase') {
    if (!isTRUE(all.equal(log, ''))) message('do you really want a logarithmic axis: log = ', log)
    for (i in (1:k)) {
      x[[i]][x[[i]] > 0.5] = x[[i]][x[[i]] > 0.5] - 1
      y[[i]] = Arg(y[[i]])/(2*pi)
      o = order(x[[i]])
      x[[i]] = sampling_rate*x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if ( (!is.null(main)) && any(is.na(main)) ) main = 'frequency response'
    if ( (!is.null(xlab)) && any(is.na(xlab)) ) {
      if (sampling_rate == 1) {
        xlab = 'frequency'
      } else {
        xlab = paste('frequency (sampling rate = ', sampling_rate, unit, ')', sep = '')
      }
    }
    if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = 'phase'
    if (any(is.na(xlim))) {
      xlim = c(0, 0.5)*sampling_rate
    }
    if (any(is.na(ylim))) {
      ylim = 'global'
    }
  }

  if (which == 'nyquist') {
    if (!isTRUE(all.equal(log, ''))) message('do you really want a logarithmic axis: log = ', log)
    for (i in (1:k)) {
      d = dim(y[[i]])
      # "close" the curves
      x[[i]] = Re(y[[i]])[,,c(1:d[3],1), drop = FALSE]
      y[[i]] = Im(y[[i]])[,,c(1:d[3],1), drop = FALSE]
    }
    if ( (!is.null(main)) && any(is.na(main)) ) main = 'Nyquist plot'
    if ( (!is.null(xlab)) && any(is.na(xlab)) )  xlab = 'real part'
    if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = 'imaginary part'
    if (any(is.na(xlim))) {
      xlim = 'subfig'
    }
    if (any(is.na(ylim))) {
      ylim = 'subfig'
    }
  }

  if ( (!is.null(subfigure_main)) && any(is.na(subfigure_main)) ) {
    if ((m0*n0) > 1) {
      subfigure_main = '(i_,j_)-th entry'
    } else {
      subfigure_main = NULL
    }
  }

  subfigure = plot_3D(x, y,
                      xlim = xlim, ylim = ylim, log = log,
                      main = main, xlab = xlab, ylab = ylab,
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      cex.points = 1, bg.points = 'black',
                      legend = legend, legend_args = legend_args)

  return(invisible(subfigure))
}


# internal helper: Spectral density to coherence
spd2coherence = function(s) {
  m = dim(s)[1]
  n.f = dim(s)[3]
  if (m == 1) {
    c = Mod(s)
    c = c/max(c)
    return(c)
  }

  c = array(0, dim = c(m,m,n.f)) # make sure that the result is a 'real' matrix!!!
  cmax = 0
  for (i in (1:m)) {
    for (j in (1:m)) {
      if (i == j) {
        # absolute vaue (however note, diagonal elements should be nonnegative)
        cmax = max(cmax, Mod(s[i,i,]))
        c[i,i,] = Mod(s[i,i,])
      }
      if (j>i) {
        # coherence
        c[i,j,] = (Mod(s[i,j,])*Mod(s[j,i,])) / (Mod(s[i,i,])*Mod(s[j,j,]))
      }
      if (j<i) {
        # phase
        c[i,j,] = Arg(s[i,j,])/(2*pi)
        # map (-0.5, 0.5) to (0, 1)
        c[i,j, c[i,j,] < 0] = c[i,j, c[i,j,] < 0] + 1
      }
    }
  }
  for (i in (1:m)) {
    # scale diagonal elements
    c[i,i,] = c[i,i,] / cmax
  }
  return(c = c)
}

#' @rdname plot
#' @export
plot.spectrald = function(x, x_list = NULL, sampling_rate = 1, unit = '',
                          which = c('modulus','phase','coherence'),
                          xlim = c(0, 0.5)*sampling_rate, ylim = 'row', log = '',
                          main = NA, xlab = NA, ylab = NA,
                          subfigure_main = NA, parse_subfigure_main = FALSE,
                          style = c('gray', 'bw', 'bw2', 'colored'),
                          col = NA, type = 'l', lty = 'solid', lwd = 1,
                          pch = 16, cex.points = 1, bg.points = 'black',
                          legend = NULL, legend_args = NA, ...) {
  style = match.arg(style)
  which = match.arg(which)

  spd = unclass(x$spd)
  attr(spd, 'z') = NULL
  m0 = dim(spd)[1]
  n.f = dim(spd)[3]
  if ((m0*n.f) == 0) stop('empty spectral density')

  x = list( (0:(n.f-1))/n.f )
  y = list(spd)

  if ((!is.null(x_list)) && (is.list(x_list))) {
    for (i in (1:length(x_list))) {
      if (!inherits(x_list[[i]],'spectrald')) stop('argument "x_list" must be a list of "spectrald" objects')

      spd = unclass(x_list[[i]]$spd)
      attr(spd, 'z') = NULL
      m = dim(spd)[1]
      n.f = dim(spd)[3]
      if (m != m0)  stop('spectrald objects are not compatible')
      if (n.f == 0) stop('empty spectral density')

      x = c(x, list( (0:(n.f-1))/n.f ) )
      y = c(y, list(spd))
    }
  }
  k = length(x)

  if ( (!is.null(main)) && (!is.expression(main)) && any(is.na(main)) ) main = 'spectral density'
  if ( (!is.null(xlab)) && (!is.expression(main))  && any(is.na(xlab)) ) {
    if (sampling_rate == 1) {
      xlab = 'frequency'
    } else {
      xlab = paste('frequency (sampling rate = ', sampling_rate, unit, ')', sep = '')
    }
  }
  if (any(is.na(xlim))) {
    xlim = c(0, 0.5)*sampling_rate
  }

  if (which == 'modulus') {
    for (i in (1:k)) {
      x[[i]][x[[i]] > 0.5] = x[[i]][x[[i]] > 0.5] - 1
      y[[i]] = abs(y[[i]])
      o = order(x[[i]])
      x[[i]] = sampling_rate * x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = 'modulus'
    if (any(is.na(ylim))) {
      ylim = 'row'
    }
  }

  if (which == 'phase') {
    for (i in (1:k)) {
      x[[i]][x[[i]] > 0.5] = x[[i]][x[[i]] > 0.5] - 1
      y[[i]] = Arg(y[[i]])/(2*pi)
      o = order(x[[i]])
      x[[i]] = sampling_rate*x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
    if ( (!is.null(ylab)) && any(is.na(ylab)) ) ylab = 'phase'
    if (any(is.na(ylim))) {
      ylim = 'row'
    }
  }

  if (which == 'coherence') {
    for (i in (1:k)) {
      x[[i]][x[[i]] > 0.5] = x[[i]][x[[i]] > 0.5] - 1
      y[[i]] = spd2coherence(y[[i]])
      o = order(x[[i]])
      x[[i]] = sampling_rate*x[[i]][o]
      y[[i]] = y[[i]][ , , o, drop = FALSE]
    }
  }

  if ((!is.null(subfigure_main)) && is.na(subfigure_main)) {
    if (m0 > 1) {
      if (which == 'coherence') {
          labels = matrix(character(m0^2), nrow = m0, ncol = m0)
          for (i in (1:m0)) {
            for (j in (1:m0)) {
              if (i == j) labels[i,j] = paste('kappa[',i,'*',j,']',sep = '')
              if (i < j)  labels[i,j] = paste('R[',i,'*',j,']^2',sep = '')
              if (i > j)  labels[i,j] = paste('Phi[',i,'*',j,']/2*pi',sep = '')
            }
          }
          subfigure_main = labels
          parse_subfigure_main = TRUE
      } else {
        subfigure_main = '(i_,j_)-th entry'
      }
    } else {
      subfigure_main = NULL
    }
  }

  # cat('plot_3D \n')
  subfigure = plot_3D(x, y,
                      xlim = xlim, ylim = ylim, log = log,
                      main = main, xlab = xlab, ylab = ylab,
                      subfigure_main = subfigure_main, parse_subfigure_main = parse_subfigure_main,
                      style = style,
                      col = col, type = type, lty = lty, lwd = lwd, pch = pch,
                      cex.points = cex.points, bg.points = bg.points,
                      legend = legend, legend_args = legend_args)
  return(invisible(subfigure))
}


#' Plot Forecast Error Variance Decomposition
#'
#' @param x [fevardec()] object.
#' @param main (character string) overall title for the plot. `main=NULL` omits the title
#'             and `main = NA` sets a default title.
#' @param xlab (character string) title for the x-axis. `xlab=NULL` omits the title
#'             and `xlab = NA` sets a default x-axis title.
#' @param col  (m)-dimensional vector of colors. If `NA` then a default colormap is
#'             chosen.
#' @param y_names optional (m)-dimensional character vector with names for the components
#'                of the time-series/process.
#' @param u_names optional (m)-dimensional character vector with names for the
#'                orthogonalized shocks.
#' @param parse_names boolean. If `TRUE` then the series- and shock- names
#'            are parsed to [expression()] before plotting. See also
#'            [grDevices::plotmath()] on the usage of expressions for
#'            plot annotations.
#' @param ...  not used.
#'
#' @return
#' This plot routine returns (invisibly) a function, subfig say, which may be used to add
#' additional graphic elements to the subfigures. The call opar = subfig(i) creates a
#' new (sub) plot at the (i)-th position with suitable margins and axis limits.
#' See the example below.
#'
#' @export
#'
#' @examples
#' # set seed for reproducible results
#' set.seed(1995)
#'
#' model = test_stspmod(dim = c(2,2), s = 3, bpoles = 1, bzeroes = 1)
#' model$names = c('x[t]', 'y[t]')
#'
#' # impulse response
#' irf = impresp(model, lag.max = 11, H = 'eigen')
#'
#' # forecast error variance decomposition
#' fevd = fevardec(irf)
#' # plot it
#' subfig = plot(fevd, col = c('lightgray','darkgray'),
#'               u_names = c('epsilon[t]', 'eta[t]'), parse_names = TRUE)
#'
#' opar = subfig(1)
#' graphics::text(x = 1, y = 0.5, 'EXAMPLE PLOT', col = 'blue', adj = c(0, 0.5))
#' graphics::par(opar)
#'
#' # reset seed
#' set.seed(NULL)

plot.fevardec = function(x, main = NA, xlab = NA, col = NA,
                         y_names = x$names, u_names = y_names,
                         parse_names = FALSE, ...) {
  vd = unclass(x)$vd
  m = dim(vd)[1]

  if (m <= 1) {
    stop('forecast error variance decomposition for a scalar process does not make sense')
  }
  if (dim(vd)[2] != m) stop('invalid "fevardec" object')

  h.max = dim(vd)[3]
  if ((m*h.max) ==0) stop('"fevardec" object is empty')

  # compute cummulative sums
  vd = apply(vd, MARGIN = c(1,3), FUN = function(x) c(0, cumsum(x)))

  # x/y limits
  xlim = c(0.45, h.max+0.55)
  ylim = c(-0.05, 1.05)

  y_names = force(y_names)
  if ( is.null(y_names) || any(is.na(y_names)) ) {
    y_names = paste('y[', 1:m, ']', sep = '')
  }
  if ( (!is.character(y_names)) || (!is.vector(y_names)) || (length(y_names) != m) ) {
    stop('illegal parameter "y_names"')
  }
  u_names = force(u_names)
  if ( is.null(u_names) || any(is.na(u_names)) ) {
    u_names = paste('u[', 1:m, ']', sep = '')
  }
  if ( (!is.character(u_names)) || (!is.vector(u_names)) || (length(u_names) != m) ) {
    stop('illegal parameter "u_names"')
  }
  u_legend = u_names
  if (parse_names) {
    u_legend = try(parse(text = u_names), silent = TRUE)
    if (inherits(u_legend, 'try-error')) u_legend = u_names
  }

  if ( (!is.null(main)) && (is.na(main)) ) main = 'forecast error variance decomposition'
  if ( (!is.null(xlab)) && (is.na(xlab)) ) xlab = 'forecast horizon'

  # get style parameters
  style = style_parameters('bw2')
  style$col.box = grDevices::gray(0.8)
  style$col.grid = grDevices::gray(0.8)
  style$lty.grid = 'solid'
  style$ldw.grid = 1
  style$bg.labels = grDevices::gray(0.8)
  style$border.labels = grDevices::gray(0.8)

  # colors
  # ()||() and ()&&() conditions must have length 1
  if (is.null(col) || all(is.na(col))) col = default_colmap(m)
  col = as.vector(col)
  col = col[(0:(m-1)) %% length(col) + 1] # recycle col if necessary

  legend_args = list(legend = u_legend, title = 'shocks',
                     fill = col, border = col, bty = 'n')

  tick_labels = array(FALSE, dim = c(m,1,4))
  tick_labels[m,1,1] = TRUE
  tick_labels[ ,1,2] = TRUE
  axes = tick_labels
  axes[,1,1] = TRUE
  titles = array(NA_character_, dim = c(m,1,4))
  titles[,1,4] = y_names

  # compute margins
  margins = list(bottom = 0.2, left = 0.2, top = 0.2, right = 0.2)
  for (i in (1:4)) {
    margins[[i]] = margins[[i]] +
      pmax(1.2 * apply(matrix(tick_labels[,,i], nrow = m, ncol = 1), MARGIN = (i-1) %% 2 + 1, FUN = any),
           0.2 * apply(matrix(axes[,,i], nrow = m, ncol = 1), MARGIN = (i-1) %% 2 + 1, FUN = any),
           1 * apply(matrix(!is.na(titles[,,i]), nrow = m, ncol = 1), MARGIN = (i-1) %% 2 + 1, FUN = any))
  }

  # set default graphic parameters
  opar = set_default_par(m,1)

  # start the plot
  start_plot(xlab = xlab, main = main, legend_args = legend_args)

  # compute the corresponding subfigures layout
  subfigure_layout = get_subfigure_layout(margins)

  subfigure_fun = function(i = 1) {
    opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,1,],
                         xaxs= 'i', yaxs = 'i',
                         fig = subfigure_layout$fig[i,1,], new = TRUE)
    # starting a new plot with plot.window() does not work ????
    graphics::plot(xlim, ylim, type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                   xlab = NA, ylab = NA, main = NA, sub = NA)
    # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
    # box(which = 'figure', col = 'red')
    return(invisible(opar))
  }

  for (i in (1:m)) {
    subfigure_fun(i)

    plot_axes(axes = axes[i,1,], tick_labels = tick_labels[i,1,], x_date = 0,
              titles = titles[i,1,], style = style, parse_titles = parse_names)

    for (j in (1:m)) {
      graphics::rect(xleft = (1:h.max) - 0.425, xright = (1:h.max) + 0.425,
                     ybottom = vd[j,i, ], ytop = vd[j+1,i,],
                     col = col[j], border = col[j])
    }
  }

  # reset graphic parameters
  graphics::par(opar)
  return(invisible(subfigure_fun))
}
