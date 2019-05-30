
#' @export
plot.CIF <- function(x, xlab ="TIME", ylab ="CUMULATIVE INCIDENCE FUNCTION",
                     xlim = NULL, ylim = NULL, main = NULL,
                     type = "l", lty = 1, lwd = 1,  col = NULL, legend = FALSE){

  if (is.null(x)){
    cat("object must be CIF.")
    return(NULL)
  }

  n <- length(x)

  max.x <- max(sapply(1:n, function(i) x[[i]][,1]))
  min.x <- min(sapply(1:n, function(i) x[[i]][,1]))
  if(is.null(xlim)) xlim <- c(min.x, max.x)

  if(is.null(col))
    col <-  grDevices::rainbow(length(n))

  if(n == 1){
    graphics::matplot(x = x[[1]][,1], y = x[[1]][,2], type = type, lty = lty,  lwd = lwd, col = col, xlab = xlab, ylab = ylab)

  }else{
    graphics::matplot(x = x[[1]][,1], y = x[[1]][,2], type = type, lty = lty,  lwd = lwd, col = col[1], xlab = xlab, ylab = ylab)

    for(i in 2:n){
      graphics::matlines(x = x[[i]][,1], y = x[[i]][,2], type = type, lty = lty,  lwd = lwd, col = col[i])

    }
  }
  legend("topleft", legend = 1:n, col = rainbow(n), lty = 1,
         title = "Individual", bty="n")

}
