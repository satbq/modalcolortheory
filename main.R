globaledo <- 12
tiny <- 1e-12

convert <- function(x,edo1,edo2) x*(edo2/edo1)
tn <- function(set, n,edo=globaledo) ((set%%edo) + (n%%edo)) %% edo
tni <- function(set, n, edo=globaledo) sort(((n%%edo) - (set%%edo)) %% edo )
startzero <- function(set, edo=globaledo) tn(set, -set[1], edo)

rotate <- function(x, n=1) {
  len <- length(x)
  n <- n %% length(x)
  return( c( tail(x,len-n), head(x,n) ))
}

rotatewrap <- function(n,x) rotate(x,n)

modecompare <- function(set, ref) sum(unique(sign(set - ref)))
# Using voice-leading brightness, modecompare returns 1 if set is brighter than ref(erence),
# -1 if set is darker than ref, and 0 if the sets are "tied" because they are identical or incomparable.

findmodes <- function(set, edo=globaledo) {
  res <- sapply(0:(length(set)-1), rotatewrap, x=set)
  res <- apply(res, 2, startzero, edo)
  return(res)
}

brightnessComps <- function(set, edo=globaledo) {
  modes <- findmodes(set, edo)
  modes <- split(modes,col(modes))
  res <- outer(modes,modes,Vectorize(modecompare))
  return(res)
}