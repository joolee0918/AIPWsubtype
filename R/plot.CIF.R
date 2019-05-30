#' Plot the cumulative incidence function
#'
#' Plotting the cumulative incidence function for an object class \code{cif}.
#'
#' @param x an object class \code{cif}.
#' @param xlab a x-axis label. Default is \code{"TIME"}.
#' @param ylab a y-axis label. Default is \code{"CUMULATIVE INCIDENCE FUNCTION"}.
#' @param main a main label.
#' @param xlim boundary for x values.
#' @param ylim boundary for y values.
#' @param type a type of plots from \code{\link[graphics]{matplot}}.
#' @param lty,lwd,lend,pch,col see \code{\link[graphics]{matplot}}.
#' @param legend a logical value: if \code{TRUE}, legend is created at the top-left corner of the plot.
#'
#'
#'
#' @importFrom graphics matplot matlines
#' @importFrom grDevices rainbow
#' @export
plot.cif <- function(x, xlab ="TIME", ylab ="CUMULATIVE INCIDENCE FUNCTION", main = NULL,
                     xlim = NULL, ylim = NULL,
                     type = "l", lty = 1, lwd = 1,  lend = par("lend"), pch = NULL, col = NULL, legend = FALSE, ...){

  if (is.null(x)){
    cat("object must be CIF.")
    return(NULL)
  }

  n <- length(x)

  max.x <- max(unlist(lapply(1:n, function(i) x[[i]]$time)))
  min.x <- min(unlist(lapply(1:n, function(i) x[[i]]$time)))

  max.y <- max(unlist(lapply(1:n, function(i) x[[i]]$cif)))
  min.y <- min(unlist(lapply(1:n, function(i) x[[i]]$cif)))

  if(is.null(xlim)) xlim <- c(min.x, max.x)
  if(is.null(ylim)) ylim <- c(min.y, max.y)

  if(is.null(col)) col <-  grDevices::rainbow(n)

  if(n == 1){
    graphics::matplot(x = x[[1]]$time, y = x[[1]]$cif, type = type, lty = lty,  lwd = lwd, lend = lend, pch = pch, col = col, xlab = xlab, ylab = ylab, ...)

  }else{
    graphics::matplot(x = x[[1]]$time, y = x[[1]]$cif, type = type, lty = lty,  lwd = lwd, lend = lend, pch = pch, col = col[1], xlab = xlab, ylab = ylab, xlim = xlim, ylim = ylim, ...)

    for(i in 2:n){
      graphics::matlines(x = x[[i]]$time, y = x[[i]]$cif, type = type, lty = lty,  lwd = lwd, lend = lend, pch = pch, col = col[i], ...)

    }
  }
  legend("topleft", legend = 1:n, col = col, lty = 1,lend = lend, pch = pch,
         title = "Individual", bty="n")

}
