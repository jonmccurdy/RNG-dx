#' DX_init
#'
#' @param K number of random numbers wanted
#' @param S number of random numbers wanted
#'
#' @return
#' @export
#'
#' @examples
#' dx_init()
#' 

dx_init <- function(K=47, S=1) {
  RNGkind("user")
  .Random.seed[2] <- K
  .Random.seed[3] <- S
  .Random.seed <<- as.integer(.Random.seed)
  invisible(runif(1))
  set.seed(.Random.seed[4])
}