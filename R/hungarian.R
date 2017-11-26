#dyn.load("R/hungarian.so")

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
    
    x<-floor(t(x) * 100)
    storage.mode(x) <- "double"
    out <- .C("hungarian", n=as.integer(nc), x=x, p=integer(nc))$p + 1
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