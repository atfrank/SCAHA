dyn.load("hungarian.so")

hungarian <-
  function(x, maximum = FALSE)
  {
    if(!is.matrix(x) || any(x < 0))
      stop("x must be a matrix with nonnegative entries.")
    
    nr <- nrow(x)
    nc <- ncol(x)
    if(nr > nc)
      stop("x must not have more rows than columns.")
    if(nc > nr)
      x <- rbind(x, matrix(2 * sum(x), nc - nr, nc))
    
    if(maximum) x <- max(x) - x
    
    # multiply all values by 1000 (3 decimal places accuracy)
    out <- .C("hungarian", n=as.integer(nc), x=as.integer(round(c(t(x)) * 1000)), p=integer(nc))$p + 1
    out <- out[seq_len(nr)]
    class(out) <- "hungarian"
    out
  }

print.hungarian <-
  function(x, ...)
  {
    writeLines(c("Optimal assignment:",
                 gsub("x", " ",
                      strwrap(paste(seq_along(x), x,
                                    sep = "x=>x", collapse = ", ")))))
    invisible(x)
  }