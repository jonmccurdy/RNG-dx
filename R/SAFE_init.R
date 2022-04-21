#' SAFE_init
#'
#'
#' @return
#' @export
#'
#' @examples
#' SAFE_init()
#' 

SAFE_init <- function(K_X=4, K_S=5) {
  if((RNGkind()!="user-supplied")[1]) RNGkind("user")
  .Random.seed[2] <- 47
  .Random.seed[3] <- 13
  .Random.seed[4] <- K_X
  .Random.seed[5] <- K_S
  .Random.seed <<- as.integer(.Random.seed)
  invisible(runif(1))
  set.seed(.Random.seed[6])
}