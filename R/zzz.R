.onLoad <- function(libname, pkgname) {
  nimble_ifelse <- nimble::nimbleFunction(
    run = function(cond = logical(0),
                   yes  = double(0),
                   no   = double(0)) {
      returnType(double(0))
      if (cond) {
        return(yes)
      } else {
        return(no)
      }
    }
  )
  
  assign(".nimble_ifelse", nimble_ifelse, envir = globalenv())
}
