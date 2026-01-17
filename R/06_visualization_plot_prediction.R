# plot_prediction.R #############################################################

#' Plot Forecasts
#'
#' The function `plot_prediction` generates some standard plots of forecasts and forecast errors.
#'
#' @param pred a list with the true data and the forecasts, as produced by [predict()].
#'        One may also add a slot "`date`" with a numeric vector of indices or a vector of type
#'        `Date` or `POSIXct` which contains date/time values.
#' @param which (character string) selects the type the plot.
#' @param qu (numeric scalar or vector) determines the width of the plotted confidence intervalls. If an entry
#'        is `NA` or equal to zero then no confidence band is plotted.
#' @param col,lty optional (vectors of) colors and line styles.
#' @param style character string determines the general style of the plot (background color, grid style,
#'              axis and axis-labels colors, ...). See also [style_parameters()].
#' @param parse_names parse series names and predictor names to [expression()].
#'        See [grDevices::plotmath()].
#' @param plot (boolean) produce a plot or just return a "closure" which then produces the plot.
#' @param ... not used
#'
#' @return If `plot=TRUE` then `plot_prediction` returns (invisibly) a function,
#'         `subfig(i = 1)` say, which may be used to add additional graphic elements to the subfigures.
#'         The call `opar = subfig(i)` creates a new (sub) plot at the (i)-th position with suitable
#'         margins and axis limits. The output `opar` contains the original graphics parameters,
#'         see [graphics::par()].
#' \cr
#' If `plot=FALSE` then a function, `plotfun(xlim = NULL)` say, is returned which produces
#' the desired plot. The optional parameter `xlim = c(x1,x2)` may be used to zoom into
#' a certain time range. The function `plotfun` returns a  function/closure to add
#' further graphical elements to the plot as described above.
#' \cr
#' See also the examples below.
#'
#' @export
#'
#' @seealso The utility [zoom_plot()] may be used to interactivly zoom in and scroll
#'          such a plot (provided that the \pkg{shiny} package is installed).
#'
#' @examples
#' # set seed for random number generation, to get reproducable results
#' set.seed(1609)
#'
#' # generate a random state space model with three outputs and 4 states
#' model = test_stspmod(dim = c(3,3), s = 4, bpoles = 1, bzeroes = 1)
#'
#' # create a vector "date" with date/time info
#' date = seq(as.POSIXct('2017-01-01'), by = 15*60, length.out = 768)
#' n.obs = sum(date < as.POSIXct('2017-01-08'))
#' n.ahead = length(date) - n.obs
#'
#' # generate random data
#' data = sim(model, n.obs = n.obs, s1 = NA)
#'
#' # compute predictions
#' pred = predict(model, data$y, h = c(1, 5), n.ahead = n.ahead)
#' # add the date/time information to the list "pred"
#' pred$date = date
#'
#' # the default "predictor names" h=1, h=2, ...
#' # don't look well, when plotted as expressions
#' dimnames(pred$yhat)[[3]] = gsub('=','==',dimnames(pred$yhat)[[3]])
#'
#' # generate some plots ####################
#'
#' # a simple/compressed plot of the data
#' p.y0 = plot_prediction(pred, which = 'y0', style = 'bw',
#'                        parse_names = TRUE, plot = FALSE)
#' # p.y0()
#'
#' # a simple/compressed plot of the prediction errors
#' plot_prediction(pred, which = 'u0', parse_names = TRUE)
#'
#' # plot of the prediction errors (with 95% confidence intervalls)
#' # plot_prediction(pred, which = 'error', qu = c(2,2,2),
#' #                 parse_names = TRUE)
#'
#' # plot of the true vales and the predicted values (+ 50% confidence region
#' # for the 1-step ahead prediction and the "out of sample" predictions)
#' p.y = plot_prediction(pred, qu = c(qnorm(0.75), NA, qnorm(0.75)),
#'                       parse_names = TRUE, plot = FALSE)
#' # subfig = p.y(xlim = date[c(n.obs-20, n.obs+20)])
#' # opar = subfig(1)
#' # abline(v = mean(as.numeric(date[c(n.obs, n.obs+1)])), col = 'red')
#' # mtext(paste(' example plot:', date()), side = 1, outer = TRUE,
#' #       cex = 0.5, col = 'gray', adj = 0)
#' # graphics::par(opar) # reset the graphical parameters
#'
#' # CUSUM plot of the prediction errors
#' # plot_prediction(pred, which = 'cusum',
#' #                 style = 'gray', parse_names = TRUE)
#'
#' # CUSUM2 plot of the prediction errors
#' # plot_prediction(pred, which = 'cusum2', parse_names = TRUE)
#'
#' set.seed(NULL) # reset seed
#'
#' \dontrun{
#' # open a 'shiny-App' window, where we can zoom
#' # into the plot with the prediction(s)
#' require(shiny)
#' zoom_plot(p.y, p.y0, 'Test zoom & scroll')
#' }

plot_prediction = function(pred, which = c('prediction', 'error', 'cusum', 'cusum2', 'y0', 'u0'),
                           qu = stats::qnorm(0.95),
                           col = NULL, lty = NULL, style = c('gray','colored','bw','bw2'),
                           parse_names = FALSE, plot = TRUE, ...) {

  which = match.arg(which)
  style = match.arg(style)
  # style parameters
  style = style_parameters(style)

  y = pred$y
  m = ncol(y)
  n.obs = nrow(y)

  y_names = colnames(pred$y)
  if (is.null(y_names)) y_names = paste('y[', 1:m, ']', sep = '')

  n.predictors = dim(pred$yhat)[3]
  if ( (!is.null(dimnames(pred$yhat))) && (!is.null(dimnames(pred$yhat)[[3]])) ) {
    p_names = dimnames(pred$yhat)[[3]]
  } else {
    p_names = paste('forecast', 1:n.predictors, sep = ifelse(parse_names, '*', ' '))
  }
  y_legend = y_names
  p_legend = p_names
  if (parse_names) {
    y_legend = try(parse(text = y_names), silent = TRUE)
    if (inherits(y_legend, 'try-error')) y_legend = y_names
    p_legend = try(parse(text = p_names), silent = TRUE)
    if (inherits(p_legend, 'try-error')) p_legend = p_names
  }

  n.ahead = nrow(pred$yhat.ahead)

  date = pred$date
  if (is.null(date)) {
    date = 1:(n.obs + n.ahead)
  }
  if (length(date)!=(n.obs+n.ahead)) {
    stop('"date" is not compatible')
  }
  x = as.numeric(date)
  x_date = date[1] - as.numeric(date[1])
  x_range = range(x)

  if (which == 'y0') {  # start of 'y0' ###################
    x = x[1:n.obs]

    y_range = apply(y, MARGIN = 1, FUN = range, na.rm = TRUE)  # 2 x n.obs

    # colors and line type
    if (is.null(col)) col = default_colmap(m)
    col = as.vector(col)
    col = col[(0:(m-1)) %% length(col) + 1] # recycle if necessary
    if (is.null(lty)) lty = 'solid'
    lty = as.vector(lty)[1]

    # legend
    if (m > 1) {
      legend_args = list(legend = y_legend, title = 'series',
                         fill = col, border = col, bty = 'n')
    } else {
      legend_args = NULL
    }

    tick_labels = c(TRUE, TRUE, FALSE, FALSE)
    axes = tick_labels
    titles = rep(NA_character_, 4)

    # compute margins
    margins = list(bottom = 1.2, left = 1.2, top = 0.2, right = 0.2)

    plotfun = function(xlim = NULL) {

      if (!is.null(xlim)) {
        xlim = as.numeric(xlim)
        i1 = which.min(abs(x - xlim[1]))
        i2 = which.min(abs(x - xlim[2]))
        # cat('xlim', xlim, dim(y_range), i1, i2,'\n')
        ylim = range(y_range[,i1:i2], na.rm = TRUE)  # 2
      } else {
        xlim = x_range
        ylim = range(y_range, na.rm = TRUE)  # 2
      }

      # set default graphic parameters
      opar = set_default_par(1,1)

      # start the plot
      start_plot(xlab = 'time (index)', main = 'time series plot', legend_args = legend_args)

      # compute the corresponding subfigures layout
      subfigure_layout = get_subfigure_layout(margins)
      # print(subfigure_layout)

      subfigure_fun = function() {
        # print(subfigure_layout)
        opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar,
                             fig = subfigure_layout$fig, new = TRUE)
        # starting a new plot with plot.window() does not work ????
        graphics::plot(xlim, ylim, type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                       xlab = NA, ylab = NA, main = NA, sub = NA)
        # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
        # box(which = 'figure', col = 'red')
        return(invisible(opar))
      }

      subfigure_fun()

      plot_axes(axes = axes, tick_labels = tick_labels, x_date = x_date,
                titles = titles, style = style, parse_titles = parse_names)

      for (j in (1:m)) {
        graphics::lines(x, y[,j], type = 'l', col = col[j], lty = lty)
      }

      # graphics::par(opar)
      return(invisible(subfigure_fun))
    } # end of plotfun

    if (plot) {
      subfigure_fun = plotfun(xlim = date[c(1,n.obs)])
      return(invisible(subfigure_fun))
    } else {
      return(invisible(plotfun))
    }
  } # end of 'y0' ###################


  if (which == 'u0') {  # start of 'u0' ###################
    x = x[1:n.obs]

    y = array(y, dim = c(n.obs, m, n.predictors)) - pred$yhat    # n.obs x m x n.predictors
    y_range = apply(y, MARGIN = 1, FUN = range, na.rm = TRUE)     # 2 x n.obs

    if ((m*n.predictors) == 1) {
      if (is.null(col)) col = 'black'
      col = as.vector(col)[1]
      if (is.null(lty)) lty = 1
      lty = as.vector(lty)[1]
      legend_args = NULL
      dim(col) = c(1,1)
      dim(lty) = c(1,1)
    }
    if ( (m==1) && (n.predictors > 1) ) {
      if (is.null(col)) col = default_colmap(n.predictors)
      col = as.vector(col)
      col = col[(0:(n.predictors-1)) %% length(col) + 1] # recycle if necessary
      if (is.null(lty)) lty = 1
      lty = as.vector(lty)[1]
      # cat(m, n.predictors, col, lty)
      legend_args = list(title = 'predictor',
                         legend = p_legend,
                         fill = col, border = col,
                         bty = 'n')
      col = matrix(col, nrow = 1, ncol = n.predictors)
      lty = matrix(lty, nrow = 1, ncol = n.predictors)
    }
    if ( (m > 1) && (n.predictors == 1) ) {
      if (is.null(col)) col = default_colmap(m)
      col = as.vector(col)
      col = col[(0:(m-1)) %% length(col) + 1] # recycle if necessary
      if (is.null(lty)) lty = 1
      lty = as.vector(lty)[1]
      legend_args = list(title = 'series',
                         legend = y_legend,
                         fill = col, border = col,
                         bty = 'n')
      col = matrix(col, nrow = m, ncol = 1)
      lty = matrix(lty, nrow = m, ncol = 1)
    }
    if ( (m > 1) && (n.predictors > 1) ) {
      if (is.null(col)) col = default_colmap(m)
      col = as.vector(col)
      col = col[(0:(m-1)) %% length(col) + 1] # recycle if necessary
      if (is.null(lty)) lty = 1:n.predictors
      lty = as.vector(lty)
      lty = lty[(0:(n.predictors-1)) %% length(lty) + 1] # recycle if necessary
      if (is.character(lty)) {
        llty = c(rep('n', m+2), lty)
      } else {
        llty = c(integer(m+2), lty)
      }
      legend_args = list(legend = c('series', y_legend, 'predictor', p_legend),
                         pch = c(NA, rep(15, m), rep(NA, n.predictors + 1)), pt.cex = 2,
                         col = c(NA, col, NA, rep('black', n.predictors)), lty = llty,
                         lwd = 2, bty = 'n')
      # print(legend_args)
      col = matrix(col, nrow = m, ncol = n.predictors)
      lty = matrix(lty, nrow = m, ncol = n.predictors, byrow = TRUE)
    }

    tick_labels = c(TRUE, TRUE, FALSE, FALSE)
    axes = tick_labels
    titles = rep(NA_character_, 4)

    # compute margins
    margins = list(bottom = 1.2, left = 1.2, top = 0.2, right = 0.2)

    plotfun = function(xlim = NULL) {

      if (!is.null(xlim)) {
        xlim = as.numeric(xlim)
        i1 = which.min(abs(x - xlim[1]))
        i2 = which.min(abs(x - xlim[2]))
        ylim = range(y_range[,i1:i2], na.rm = TRUE)  # 2
      } else {
        xlim = x_range
        ylim = range(y_range, na.rm = TRUE)  # 2
      }

      # set default graphic parameters
      opar = set_default_par(1,1)

      # start the plot
      start_plot(xlab = 'time (index)', main = 'time series plot of prediction errors', legend_args = legend_args)

      # compute the corresponding subfigures layout
      subfigure_layout = get_subfigure_layout(margins)

      subfigure_fun = function() {
        opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar,
                             fig = subfigure_layout$fig, new = TRUE)
        # starting a new plot with plot.window() does not work ????
        graphics::plot(xlim, ylim, type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                       xlab = NA, ylab = NA, main = NA, sub = NA)
        # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
        # box(which = 'figure', col = 'red')
        return(invisible(opar))
      }

      subfigure_fun()

      plot_axes(axes = axes, tick_labels = tick_labels, x_date = x_date,
                titles = titles, style = style, parse_titles = parse_names)

      for (i in (1:m)) {
        for (j in (1:n.predictors)) {
          graphics::lines(x, y[,i,j], type = 'l', col = col[i,j], lty = lty[i,j])
        }
      }

      # graphics::par(opar)
      return(invisible(subfigure_fun))
    } # end of plotfun

    if (plot) {
      subfigure_fun = plotfun(xlim = date[c(1,n.obs)])
      return(invisible(subfigure_fun))
    } else {
      return(invisible(plotfun))
    }
  } # end of 'u0' ###################


  if (which == 'error') {  # start of 'error' ###################
    x = x[1:n.obs]

    # quantiles
    qu = as.vector(qu)
    qu = qu[(0 : (n.predictors - 1)) %% length(qu) + 1] # recycle if necessary

    # standard deviations times the quantiles
    sd = matrix(qu, nrow = m, ncol = n.predictors, byrow = TRUE) *
      apply(pred$sigmahat, MARGIN = 3, FUN = function(x) {sqrt(diag(matrix(x, nrow = m, ncol = m)))})
    sd = aperm( array(sd, dim = c(m, n.predictors, n.obs)), c(3,1,2) ) # n.obs x m x n.predictors


    y = array(y, dim = c(n.obs, m, n.predictors)) - pred$yhat    # n.obs x m x n.predictors

    y_range1 = apply(sd, MARGIN = c(1,2), FUN = max, na.rm = TRUE)    # n.obs x m
    dim(y_range1) = c(1, n.obs, m)                                    # 1 x n.obs x m
    y_range2 = apply(y, MARGIN = c(1,2), FUN = range, na.rm = TRUE)   # 2 x n.obs x m
    y_range = apply(dbind(d = 1, y_range2, y_range1, -y_range1), MARGIN = c(2,3),
                    FUN = range, na.rm = TRUE)                        # 2 x n.obs x m

    # line colors, line type and fill colors
    if (is.null(col)) col = default_colmap(n.predictors)
    col = as.vector(col)
    col = col[(0:(n.predictors-1)) %% length(col) + 1] # recycle if necessary
    if (is.null(lty)) lty = 'solid'
    lty = as.vector(lty)[1]
    alpha = grDevices::col2rgb(col)
    alpha = grDevices::rgb(red = alpha[1,], green = alpha[2,], blue = alpha[3,],
                alpha = rep(round(0.25*255), n.predictors), maxColorValue = 255)

    # legend
    if (n.predictors > 1) {
      legend_args = list(legend = p_legend, title = 'predictors',
                         fill = col, border = col, bty = 'n')
    } else {
      legend_args = NULL
    }

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

    plotfun = function(xlim = NULL) {

      if (!is.null(xlim)) {
        xlim = as.numeric(xlim)
        i1 = which.min(abs(x - xlim[1]))
        i2 = which.min(abs(x - xlim[2]))
        ylim = apply(y_range[,i1:i2,,drop = FALSE], MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      } else {
        xlim = x_range
        ylim = apply(y_range, MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      }

      # set default graphic parameters
      opar = set_default_par(m,1)

      # start the plot
      start_plot(xlab = 'time (index)', main = 'time series plot of prediction errors',
                 legend_args = legend_args)

      # compute the corresponding subfigures layout
      subfigure_layout = get_subfigure_layout(margins)

      subfigure_fun = function(i = 1) {
        opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,1,],
                             fig = subfigure_layout$fig[i,1,], new = TRUE)
        # starting a new plot with plot.window() does not work ????
        graphics::plot(xlim, ylim[,i], type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                       xlab = NA, ylab = NA, main = NA, sub = NA)
        # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
        # box(which = 'figure', col = 'red')
        return(invisible(opar))
      }

      for (i in (1:m)) {
        subfigure_fun(i)

        plot_axes(axes = axes[i,1,], tick_labels = tick_labels[i,1,], x_date = x_date,
                  titles = titles[i,1,], style = style, parse_titles = parse_names)

        for (j in (1:n.predictors)) {
          graphics::polygon(c(x, rev(x)), c(-sd[,i,j], sd[,i,j]),
                  col = alpha[j], border = alpha[j])
          graphics::lines(x, y[,i,j], type = 'l', col = col[j], lty = lty)
        }
      }

      # graphics::par(opar)
      return(invisible(subfigure_fun))
    } # end of plotfun

    if (plot) {
      subfigure_fun = plotfun(xlim = date[c(1,n.obs)])
      return(invisible(subfigure_fun))
    } else {
      return(invisible(plotfun))
    }
  } # end of 'error' ###################


  if (which == 'prediction') {  # start of 'prediction' ###################

    # quantiles
    qu = as.vector(qu)
    qu = qu[(0 : n.predictors) %% length(qu) + 1] # recycle if necessary  (n.predictors + 1)

    # standard deviations times the quantiles
    sd = matrix(qu[1:n.predictors], nrow = m, ncol = n.predictors, byrow = TRUE) *
      apply(pred$sigmahat, MARGIN = 3, FUN = function(x) {sqrt(diag(matrix(x, nrow = m, ncol = m)))})
    sd = aperm( array(sd, dim = c(m, n.predictors, n.obs)), c(3,1,2) ) # n.obs x m x n.predictors

    if (n.ahead > 0) {
      sd.ahead = qu[n.predictors+1] *
        apply(pred$sigmahat.ahead, MARGIN = 3, FUN = function(x) {sqrt(diag(matrix(x, nrow = m, ncol = m)))})
      sd.ahead = t(matrix(sd.ahead, nrow = m, ncol = n.ahead))  # n.ahead x m
    } else {
      sd.ahead = matrix(0, nrow = 0, ncol = m)
    }

    yhat = pred$yhat                # n.obs x m x n.predictors
    yhat.ahead = pred$yhat.ahead    # n.ahead x m

    y_range1 = apply(y, MARGIN = c(1,2), FUN = range, na.rm = TRUE)             # 2 x n.obs x m
    y_range2 = apply(yhat + sd, MARGIN = c(1,2), FUN = range, na.rm = TRUE)     # 2 x n.obs x m
    y_range3 = apply(yhat - sd, MARGIN = c(1,2), FUN = range, na.rm = TRUE)     # 2 x n.obs x m
    y_range.in = apply(dbind(d = 1, y_range1, y_range2, y_range3),
                       MARGIN = c(2,3), FUN = range, na.rm = TRUE)              # 2 x n.obs x m
    if (n.ahead > 0) {
      y_range.ahead = array(0, dim = c(2, n.ahead, m))
      y_range.ahead[1,,] = yhat.ahead - sd.ahead
      y_range.ahead[2,,] = yhat.ahead + sd.ahead
      y_range = dbind(d = 2, y_range.in, y_range.ahead)                           # 2 x (n.obs+n.ahead) x m
    } else {
      y_range = y_range.in
    }

    # line colors, line type and fill colors
    if (is.null(col)) col = default_colmap(n.predictors + 1)
    col = as.vector(col)
    col = col[(0:n.predictors) %% length(col) + 1] # recycle if necessary
    if (is.null(lty)) lty = 'solid'
    lty = as.vector(lty)[1]
    alpha = grDevices::col2rgb(col)
    alpha = grDevices::rgb(red = alpha[1,], green = alpha[2,], blue = alpha[3,],
                alpha = rep(round(0.25*255), n.predictors + 1), maxColorValue = 255)

    # legend
    if (n.predictors > 1) {
      legend_args = list(legend = c('true', p_legend, 'ahead'),
                         title = 'predictors',
                         fill = c(grDevices::gray(0.25), col),
                         border = c(grDevices::gray(0.25), col), bty = 'n')
    } else {
      legend_args = NULL
    }

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

    plotfun = function(xlim = NULL) {

      if (!is.null(xlim)) {
        xlim = as.numeric(xlim)
        i1 = which.min(abs(x - xlim[1]))
        i2 = which.min(abs(x - xlim[2]))
        ylim = apply(y_range[,i1:i2,,drop = FALSE], MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      } else {
        xlim = x_range
        ylim = apply(y_range, MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      }

      # set default graphic parameters
      opar = set_default_par(m,1)

      # start the plot
      start_plot(xlab = 'time (index)', main = 'time series plot of predictions', legend_args = legend_args)

      # compute the corresponding subfigures layout
      subfigure_layout = get_subfigure_layout(margins)

      subfigure_fun = function(i = 1) {
        opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,1,],
                             fig = subfigure_layout$fig[i,1,], new = TRUE)
        # starting a new plot with plot.window() does not work ????
        graphics::plot(xlim, ylim[,i], type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                       xlab = NA, ylab = NA, main = NA, sub = NA)
        # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
        # box(which = 'figure', col = 'red')
        return(invisible(opar))
      }

      for (i in (1:m)) {
        subfigure_fun(i)

        plot_axes(axes = axes[i,1,], tick_labels = tick_labels[i,1,], x_date = x_date,
                  titles = titles[i,1,], style = style, parse_titles = parse_names)

        for (j in (1:n.predictors)) {
          graphics::polygon(c(x[1:n.obs], x[n.obs:1]), c(yhat[,i,j] - sd[,i,j], rev(yhat[,i,j] + sd[,i,j])),
                  col = alpha[j], border = alpha[j])
          graphics::lines(x[1:n.obs], yhat[,i,j], col = col[j], lty = lty, lwd = 1)
        }
        if (n.ahead > 0) {
          graphics::polygon(c(x[(n.obs+1):(n.obs+n.ahead)], x[(n.obs+n.ahead):(n.obs+1)]),
                  c(yhat.ahead[,i] - sd.ahead[,i], rev(yhat.ahead[,i] + sd.ahead[,i])),
                  col = alpha[n.predictors+1], border = alpha[n.predictors+1])

          graphics::lines(x[(n.obs+1):(n.obs+n.ahead)], yhat.ahead[,i],
                col = col[n.predictors+1], lty = lty, lwd = 1)
        }
        graphics::lines(x[1:n.obs], y[,i], type = 'l', col = grDevices::gray(0.25), lty = lty)
      }

      # graphics::par(opar)
      return(invisible(subfigure_fun))
    } # end of plotfun

    if (plot) {
      subfigure_fun = plotfun()
      return(invisible(subfigure_fun))
    } else {
      return(invisible(plotfun))
    }
  } # end of 'prediction' ###################


  if ((which == 'cusum') || (which == 'cusum2')) { # start of 'cusum'/'cusum2' ###########
    cum2 = (which == 'cusum2')
    x = x[1:n.obs]

    y = array(y, dim = c(n.obs, m, n.predictors)) - pred$yhat
    if (cum2) {
      y = apply(y^2, MARGIN = c(2,3), FUN = cumsum)                # n.obs x m x n.predictors

      # prediction error variances m x n.predictors
      var = matrix(apply(pred$sigmahat, MARGIN = 3,
                         FUN = function(x) { diag(matrix(x, nrow = m, ncol = m)) }),
                   nrow = m, ncol = n.predictors)
      bb = var * (n.obs / (x[n.obs] - x[1]))
      aa = -x[1] * bb
    } else {
      y = apply(y, MARGIN = c(2,3), FUN = cumsum)                  # n.obs x m x n.predictors
    }
    y_range = apply(y, MARGIN = c(1,2), FUN = range, na.rm = TRUE)  # 2 x n.obs x m

    if (cum2) {
    }

    # colors and line type
    if (is.null(col)) col = default_colmap(n.predictors)
    col = as.vector(col)
    col = col[(0:(n.predictors-1)) %% length(col) + 1] # recycle if necessary
    if (is.null(lty)) lty = 'solid'
    lty = as.vector(lty)[1]

    # legend
    if (n.predictors > 1) {
      legend_args = list(legend = p_legend, fill = col, border = col, bty = 'n', title = 'predictors')
    } else {
      legend_args = NULL
    }

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

    plotfun = function(xlim = NULL) {

      if (!is.null(xlim)) {
        xlim = as.numeric(xlim)
        i1 = which.min(abs(x - xlim[1]))
        i2 = which.min(abs(x - xlim[2]))
        ylim = apply(y_range[,i1:i2,,drop = FALSE], MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      } else {
        xlim = x_range
        ylim = apply(y_range, MARGIN = 3, FUN = range, na.rm = TRUE)  # 2 x m
      }

      # set default graphic parameters
      opar = set_default_par(m,1)

      # start the plot
      start_plot(xlab = 'time (index)',
                 main = paste(ifelse(cum2, 'CUSUM2', 'CUSUM'), 'chart of prediction errors'),
                 legend_args = legend_args)

      # compute the corresponding subfigures layout
      subfigure_layout = get_subfigure_layout(margins)

      subfigure_fun = function(i = 1) {
        opar = graphics::par(omi = subfigure_layout$omi, mar = subfigure_layout$mar[i,1,],
                             fig = subfigure_layout$fig[i,1,], new = TRUE)
        # starting a new plot with plot.window() does not work ????
        graphics::plot(xlim, ylim[,i], type = 'n', log = '', axes = FALSE, frame.plot = FALSE,
                       xlab = NA, ylab = NA, main = NA, sub = NA)
        # cat('subfig_fun:', par('oma'), '\n', mar, '\n', fig, '\n', xlim, '\n', ylim, '\n')
        # box(which = 'figure', col = 'red')
        return(invisible(opar))
      }

      for (i in (1:m)) {
        subfigure_fun(i)
        # cat(i,'\n')
        # graphics::lines(xlim, ulim[,i])
        plot_axes(axes = axes[i,1,], tick_labels = tick_labels[i,1,], x_date = x_date,
                  titles = titles[i,1,], style = style, parse_titles = parse_names)

        if (cum2) {
          for (j in (1:n.predictors)) {
            graphics::abline(a = aa[i,j], b = bb[i,j], col = col[j], lty = 'dotted')
          }
        } else {
          graphics::abline(h = 0, col = 'darkgray')
        }

        # usr = par('usr')
        # rect(usr[1], usr[3], usr[2], usr[4], col = 'orange')
        for (j in (1:n.predictors)) {
          graphics::lines(x, y[,i,j], type = 'l', col = col[j], lty = lty)
        }
      }

      # box(which = 'outer')
      # box(which = 'inner', col = 'red')

      # graphics::par(opar)
      return(invisible(subfigure_fun))
    } # end of plofun

    if (plot) {
      subfigure_fun = plotfun(xlim = date[c(1,n.obs)])
      return(invisible(subfigure_fun))
    } else {
      return(invisible(plotfun))
    }
  }  # end of 'cusum'/'cusum2' ###########

  stop('this should not happen')
}

