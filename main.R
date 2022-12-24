globaledo <- 12
tiny <- 1e-12
globalrounder <- -log(tiny)/log(10)
fortenums <- readRDS("fortenums.rds")

fpunique <- function(x, MARGIN=0, rounder=globalrounder) {
  if (MARGIN == 0) {
    return(x[!duplicated(round(x, rounder))])
  }

  if (MARGIN == 1) {
    return(x[!duplicated( round(x, rounder), MARGIN=1 ) , ] )
  }

  if (MARGIN == 2) {
    return(x[ , !duplicated( round(x, rounder), MARGIN=2 ) ] )
  }
}

convert <- function(x, edo1, edo2) x*(edo2/edo1)

tn <- function(set, n, edo=globaledo, sorted=TRUE) {
  res <- ((set%%edo) + (n%%edo)) %% edo
  if (sorted == FALSE) { return(res) }
  return(sort(res))
}

tni <- function(set, n, edo=globaledo, sorted=TRUE) {
  res <- ((n%%edo) - (set%%edo)) %% edo
  if (sorted == FALSE) { return(res) }
  return(sort(res))
}

startzero <- function(set, edo=globaledo, sorted=TRUE) tn(set, -set[1], edo, sorted)

rotate <- function(x, n=1) {
  len <- length(x)
  n <- n %% length(x)
  return( c( tail(x,len-n), head(x,n) ))
}

rotatewrap <- function(n,x) rotate(x,n)

modecompare <- function(set, ref, rounder=globalrounder) sum(unique(sign(round(set - ref, rounder))))
# Using voice-leading brightness, modecompare returns 1 if set is brighter than ref(erence),
# -1 if set is darker than ref, and 0 if the sets are "tied" because they are identical or incomparable.

# sim = scalar interval matrix
sim <- function(set, edo=globaledo) {
  res <- sapply(0:(length(set)-1), rotatewrap, x=set)
  res <- apply(res, 2, startzero, edo)
  return(res)
}

brightnessComps <- function(set, edo=globaledo, rounder=globalrounder) {
  modes <- sim(set, edo)
  modes <- split(modes,col(modes))
  res <- outer(modes,modes,Vectorize(modecompare), rounder=rounder)
  return(res)
}

eps <- function(set, edo=globaledo, rounder=globalrounder) {
  modes <- t(sim(set, edo))
  chart <- brightnessComps(set, edo, rounder)*brightnessComps(set, edo, rounder)
  diffs <- outer(rowSums(modes),rowSums(modes),'-')
  result <- chart * diffs

  return(min(result[result > 0]))
}

delta <- function(set, edo=globaledo, rounder=globalrounder) {
  modes <- t(sim(set, edo))
  chart <- brightnessComps(set, edo, rounder)*brightnessComps(set, edo, rounder)
  diag(chart) <- -1
  chart <- (chart + 1)%%2

  diffs <- outer(rowSums(modes),rowSums(modes),'-')
  result <- chart * diffs

  return(max(result))
}

ratio <- function(set, edo=globaledo, rounder=globalrounder) {
  return(delta(set, edo, rounder)/eps(set, edo, rounder))
}

sc <- function(card,num) {
  set <- fortenums[[card]][num]
  res <- strtoi(unlist(strsplit(set,split=",")))
  return(res)
}

tnprime <- function(set, edo=globaledo) {
  set <- sort(set %% edo)
  card <- length(set)
  if (card == 1) { return(0) }
  if (card == 0) { return(integer(0))}
  modes <- sim(set, edo)

  for (i in card:1) {
     if (class(modes)[1]=="numeric") {return(modes)}
     top <- min(modes[i,])
     index <- which(modes[i,] == top)
     modes <- modes[,index]
  }

  return(modes[,1])

}

setcompare <- function(x,y) {
  card <- length(x)
  if ( length(y) != card ) { print("Cardinality mismatch"); return(NA) }

  modes <- cbind(x,y)

  for (i in card:1) {
    if (class(modes)[1]=="numeric") {return(modes)}
    top <- min(modes[i,])
    index <- which(modes[i,] == top)
    modes <- modes[,index]
  }

  return(modes[,1])
}

primeform <- function(set, edo=globaledo) {
  if (length(set)==1) { return(0) }
  upset <- startzero(tnprime(set, edo))
  downset <- startzero(tnprime(tni(set, 0, edo), edo))
  winner <- setcompare(upset, downset)
  return(winner)
}

charm <- function(set, edo=globaledo) {
  return(tnprime(tni(set, 0, edo), edo))
}

#set-class complemnent
scComp <- function(set,edo=globaledo) {
  return(primeform(setdiff(0:(edo-1),set),edo))
}

fortenum <- function(set) {
  card <- length(set)
  strset <- toString(primeform(set, edo=12))
  val <- which(fortenums[[card]]==strset)
  return(paste0(card, "-", val))
}

zmate <- function(set) {
  num <- fortenum(set)
  res <- NA

  if (num == "4-15") { res <- c(4,29) }
  if (num == "4-29") { res <- c(4,15) }

  if (num == "8-15") { res <- c(8,29) }
  if (num == "8-29") { res <- c(8,15) }

  if (num == "5-12") { res <- c(5,36) }
  if (num == "5-36") { res <- c(5,12) }
  if (num == "5-38") { res <- c(5,18) }
  if (num == "5-18") { res <- c(5,38) }
  if (num == "5-37") { res <- c(5,17) }
  if (num == "5-17") { res <- c(5,37) }

  if (num == "7-12") { res <- c(7,36) }
  if (num == "7-36") { res <- c(7,12) }
  if (num == "7-38") { res <- c(7,18) }
  if (num == "7-18") { res <- c(7,38) }
  if (num == "7-37") { res <- c(7,17) }
  if (num == "7-17") { res <- c(7,37) }


  if (num == "6-36") { res <- c(6,3) }
  if (num == "6-3") { res <- c(6,36) }
  if (num == "6-37") { res <- c(6,4) }
  if (num == "6-4") { res <- c(6,37) }
  if (num == "6-40") { res <- c(6,11) }
  if (num == "6-11") { res <- c(6,40) }
  if (num == "6-41") { res <- c(6,12) }
  if (num == "6-12") { res <- c(6,41) }
  if (num == "6-13") { res <- c(6,42) }
  if (num == "6-42") { res <- c(6,13) }
  if (num == "6-38") { res <- c(6,6) }
  if (num == "6-6") { res <- c(6,38) }
  if (num == "6-46") { res <- c(6,24) }
  if (num == "6-24") { res <- c(6,46) }
  if (num == "6-17") { res <- c(6,43) }
  if (num == "6-43") { res <- c(6,17) }
  if (num == "6-47") { res <- c(6,25) }
  if (num == "6-25") { res <- c(6,47) }
  if (num == "6-19") { res <- c(6,44) }
  if (num == "6-44") { res <- c(6,19) }
  if (num == "6-48") { res <- c(6,26) }
  if (num == "6-26") { res <- c(6,48) }
  if (num == "6-10") { res <- c(6,39) }
  if (num == "6-39") { res <- c(6,10) }
  if (num == "6-49") { res <- c(6,28) }
  if (num == "6-28") { res <- c(6,49) }
  if (num == "6-29") { res <- c(6,50) }
  if (num == "6-50") { res <- c(6,29) }
  if (num == "6-45") { res <- c(6,23) }
  if (num == "6-23") { res <- c(6,45) }

  if (is.na(res[1])) { return(res) }

  return(sc(res[1],res[2]))
}

ivec <- function(set, edo=globaledo) {
  set <- set%%edo
  set <- unique(set)
  vec <- rep(edo+1,edo/2)
  ivs <- outer(set, set, "-")
  ivs2 <- (edo - ivs)
  lowers <- ivs
  lowers[which(ivs > ivs2)] <- ivs2[which(ivs > ivs2)]
  nonzero <- lowers[lowers > 0]

  for (i in 1:(edo/2)) {
    vec[i] <- sum(nonzero == i)
  }

  return(vec)
}

# "equal division of the octave origin"
edoo <- function(card, edo=globaledo) {
  return( (0:(card-1))*(edo/card) )
}

evenness <- function(set, edo=globaledo) {
  card <- length(set)
  edoozero <- edoo(card, edo) - (sum(edoo(card, edo))/card)
  setzero <- set - (sum(set)/card)
  return(sqrt(sum((setzero-edoozero)^2)))
}

saturate <- function(r, set, edo=globaledo) {
  card <- length(set)
  origin <- edoo(card, edo)
  vec <- origin + r*(set-origin)

  return(vec)
}

asword <- function(set, edo=globaledo, rounder=globalrounder) {
  card <- length(set)
  result <- rep(0,card)

  setsteps <- sim(set, edo)[2, ]
  stepvals <- sort(fpunique(setsteps, 0, rounder))

  for (i in 1:card) {
    result[which(abs(setsteps - stepvals[i]) < 10^(-1*rounder) )] <- i
  }

  return(result)
}

intervalspectrum <- function(set, edo=globaledo, rounder=globalrounder) {
  modes <- sim(set, edo)
  uniques <- apply(modes,1,fpunique, rounder=rounder)
  return(uniques)
}

spectrumcount <- function(set, edo=globaledo, rounder=globalrounder) sapply(intervalspectrum(set,edo,rounder),length)

diatonicsubsets <- function(subsetdegrees,set,uniq=TRUE,edo=globaledo,rounder=globalrounder) {
  modes <- sim(set,edo=edo)
  subsetdegrees <- subsetdegrees + 1 #Because in music theory these are 0-indexed, but vectors are 1-indexed in R
  res <- modes[subsetdegrees,]

  if (uniq == TRUE) {
  res <- fpunique(res, MARGIN=2, rounder=rounder)
  }

  return(res)
}

subsetspectrum <- function(set,subsetcard,simplify=TRUE,mode="tn",edo=globaledo,rounder=globalrounder) { #Returns a list
  # edosave <- edo
  card <- length(set)
  comb <- combn(card-1,subsetcard-1)
  comb <- rbind(rep(0,choose(card-1,subsetcard-1)),comb)

  if (mode=="tn") { use <- tnprime }
  if (mode=="tni") { use <- primeform }

  if (simplify == TRUE) {
    # edo <<- card
    comb <- fpunique(apply(comb,2,use,edo=card),MARGIN=2,rounder=rounder)
    # edo <<- edosave
  }

  res <- apply(comb,2,diatonicsubsets,set=set,edo=edo,rounder=rounder)

  if ("matrix" %in% class(res)) {	#This is necessary because if all subsets have same variety, apply() returns a matrix, not a list
    res <- as.list(as.data.frame(res))

    for (i in 1:length(res) ) {
       res[[i]] <- matrix(res[[i]],nrow=subsetcard,ncol=(length(res[[i]])/subsetcard))
    }
  }

  names(res) <- apply(comb,2,toString)
  return(res)
}

isym <- function(set, edo=globaledo, rounder=globalrounder) {
  card <- length(set)
  setword <- asword(set, edo, rounder)
  invsetword <- rev(setword)

  for (i in 1:card) {
    invmode <- rotate(invsetword,i)
    if ( isTRUE(all.equal(setword,invmode)) ) { return(TRUE) }
  }
  return(FALSE)
}