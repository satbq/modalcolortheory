globaledo <- 12
tiny <- 1e-12
globalrounder <- -log(tiny)/log(10)
fortenums <- readRDS("fortenums.rds")
ineqmats <- readRDS("ineqmats.rds")

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


makesubmatold <- function(card) {
  submat <- matrix("",card,card)

  for (i in 0:(card-1)) {
    submat[(i+1),] <- rep(paste0("-",i),card)
  }

  for (i in 1:card) {
    prefixes <- 0:(card-1)
    prefixes <- rotate(prefixes,(i-1))

    for (j in 1:card) {
      submat[i,j] <- paste0(prefixes[j], submat[i,j])
    }
  }

  for (i in 1:card) {
    for (j in 1:card) {
      if ((i+j)>(card+1)) { submat[i,j] <- paste0("w+", submat[i,j])}
    }
  }
  return(submat)
}

#MUCH faster than the original loop-based version
makesubmat <- function(card) {
  submat <- t(outer(0:(card-1), 0:(card-1), paste, sep="-"))

  rotaterow <- function(n,mat) rotate(mat[n,], n-1)
  submat <- t(sapply(1:card, rotaterow,mat=submat))

  submat <- submat[,card:1]
  submat[lower.tri(submat)] <- mapply(paste0, "w+", submat[lower.tri(submat)])
  submat <- submat[,card:1]

  return(submat)
}

makeineqmatold <- function(card) {
  # This is fairly slow. Not intended for repeated use: look up precalculated results with getineqmat!
  ineqmat <- rep(0,(card+1))
  ineqmat <- rbind(ineqmat,ineqmat)
  submat <- t(makesubmat(card))

  for (i in 2:card) {
    for (j in 1:(card-1)) {
      posvec <- submat[i,j]

      for (z in (j+1):card) {
        negvec <- submat[i,z]
        newline <- rep(0,(card+1))
        wtest <- 1

        if (substr(posvec,1,1)=="w") {
          wtest <- wtest*(-1)
          q <- nchar(posvec)
          posvec <- substr(posvec,3,q)
        }

        if (substr(negvec,1,1)=="w") {
          wtest <- wtest*(-1)
          q2 <- nchar(negvec)
          negvec <- substr(negvec,3,q2)
        }

        q <- nchar(posvec)
        q2 <- nchar(negvec)

        if (wtest == (-1)) { newline[card+1] <- -1 }

        substrRight <- function(x, n){
          # This function is by user Andrie on Stack Overflow: https://stackoverflow.com/a/7963963
          substr(x, nchar(x)-n+1, nchar(x))
        }

        q3 <- as.integer(gregexpr(pattern="-",posvec)[[1]])
        q4 <- as.integer(gregexpr(pattern="-",negvec)[[1]])
        posvec1 <- as.numeric(substr(posvec,1,(q3-1))) + 1
        posvec2 <- as.numeric(substrRight(posvec,(q-q3))) + 1
        negvec1 <- as.numeric(substr(negvec,1,(q4-1))) + 1
        negvec2 <- as.numeric(substrRight(negvec,(q2-q4))) + 1

        newline[posvec1] <- newline[posvec1] + 1
        newline[negvec2] <- newline[negvec2] + 1
        newline[posvec2] <- newline[posvec2] - 1
        newline[negvec1] <- newline[negvec1] - 1

        negline <- newline
        negline[1:card] <- -1 * negline[1:card]

        newtest <- TRUE

        temp <- rbind(ineqmat[,1:card],newline[1:card])
        counter <- dim(temp)[1]
        temp2 <- unique(temp,MARGIN=1)
        counter2 <- dim(temp2)[1]

        if (counter != (counter2+1)) { newtest <- FALSE }

        temp <- rbind(ineqmat[,1:card],negline[1:card])
        counter <- dim(temp)[1]
        temp2 <- unique(temp,MARGIN=1)
        counter2 <- dim(temp2)[1]
        if (counter != (counter2+1)) { newtest <- FALSE }

        if ( newtest == TRUE ) { ineqmat <- rbind(ineqmat,newline) }

      }

    }
  }

  ineqmat <- ineqmat[-(1:2),]
  row.names(ineqmat) <- NULL
  return(ineqmat)
}

makeineqmat <- function(card) {
  # Creates a row for the inequality matrix, given the "roots" of the two intervals to be compared
  # (specified as zero-indexed scale degrees) and the generic size of the interval (zero-indexed).
  generateRow <- function(firstroot, secondroot, genericival) {
    row <- rep(0, card+1)
    if ((secondroot %% card) <= firstroot) { return(row) }
    row[(firstroot %% card)+1] <- row[(firstroot %% card)+1] - 1
    row[(secondroot %% card)+1] <- row[(secondroot %% card)+1] + 1
    row[((firstroot + genericival) %% card) + 1] <- row[((firstroot + genericival) %% card) + 1] + 1
    row[((secondroot + genericival) %% card)+1] <- row[((secondroot + genericival) %% card)+1] - 1

    # Last column of the inequality matrix reflects whether the intervals in the comparison wrap around the octave.
    # For instance, comparing do-re to ti-do requires adding 12 to the higher do.
    w <- ((firstroot + genericival) >= card) - ((secondroot + genericival) >= card)
    row[card+1] <- w
    return(row)
  }

  roots <- 0:(card-1)
  intervals <- 1:(card/2)

  # Create all possible combinations of roots and interval sizes.
  combinations <- expand.grid(roots, roots, intervals)
  firstroots <- combinations[,1]
  secondroots <- combinations[,2]
  genericintervals <- combinations[,3]

  # Generate rows from all the combinations; many will be redundant.
  res <- t(mapply(generateRow, firstroot=firstroots, secondroot=secondroots, genericival=genericintervals))

  # Next two lines guarantee that we'll only generate each hyperplane in one orientation.
  # (We don't want to have separate hyperplanes for "first step bigger than second step" and "first step smaller...")
  # So arbitrarily we require the first nonzero entry of a row to be negative.
  rowSign <- function(row) row * -1 * sign(row[which(row!=0)])[1]
  res <- t(apply(res, 1, rowSign))

  # Remove redundant rows.
  res <- res[!duplicated(res, MARGIN=1),]

  # First row was all 0s -- unisons are identical -- so remove it.
  res <- res[-1,]

  return(res)
}

if (FALSE) {
  # The code that generates ineqmats.rds. No need to run!
  ineqmats <- list()
  ineqmats[[1]] <- integer(0)
  for (i in 2:53) {
    time1 <- Sys.time()
    ineqmats[[i]] <- makeineqmat(i)
    time2 <- Sys.time()
    print(i)
    print(time2-time1)
  }
  saveRDS(ineqmats, "ineqmats.rds")
}

getineqmat <- function(card) {
  if (exists("ineqmats")) {
    if (card <= length(ineqmats)) {
      return(ineqmats[[card]])
    }
  }
 else {
    return(makeineqmat(card))
  }
}

signvector <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  if (is.null(ineqmat)) {
    card <- length(set)
    ineqmat <- getineqmat(card)
  }
  set <- c(set, edo)
  res <- ineqmat %*% set
  res <- sign(round(res, digits=rounder))
  return(as.vector(res))
}

whichsvzeroes <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  signvec <- signvector(set, ineqmat, edo, rounder)
  return(which(signvec == 0))
}

countsvzeroes <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  signveczeroes <- whichsvzeroes(set, ineqmat, edo, rounder)
  return(length(signveczeroes))
}

howfree <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  card <- length(set)

  if (is.null(ineqmat)) {
    ineqmat <- getineqmat(card)
  }

  zeroesflat <- ineqmat[whichsvzeroes(set, ineqmat, edo, rounder),]
  rank <- qr(zeroesflat)$rank
  freedom <- card - (1+rank)
  return(freedom)
}

# Checks relationship between two signvecs. If identical returns 0; if adjacent regions, returns 1; otherwise returns -1.
comparesignvecs <- function(signvecX, signvecY) {
  # Are the signvecs identical?
  if ( isTRUE(all.equal(signvecX,signvecY)) ) { return(0) }

  # Initial sorting based on which hyperplanes the scales lie on.
  signvec.zeroes.x <- which(signvecX == 0)
  signvec.zeroes.y <- which(signvecY == 0)

  numzeroes.x <- length(signvec.zeroes.x)
  numzeroes.y <- length(signvec.zeroes.y)

  # Distinct colors can't be adjacent if lie on same number of hyperplanes.
  if (numzeroes.x == numzeroes.y) {
    return(-1)
  }

  # Which scale lies on more hyperplanes?
  if (numzeroes.x > numzeroes.y) {
    upper <- signvecX
    upper.zeroes <- signvec.zeroes.x
    lower <- signvecY
    lower.zeroes <- signvec.zeroes.y
  }

  if (numzeroes.x < numzeroes.y) {
    upper <- signvecY
    upper.zeroes <- signvec.zeroes.y
    lower <- signvecX
    lower.zeroes <- signvec.zeroes.x
  }

  # To be adjacent, a freer region must lie on all the hyperplanes of the stricter region.
  check.zeroes <- lower.zeroes %in% upper.zeroes
  if ( prod(check.zeroes) == 0 ) {
    return(-1)
  }

  # We've weeded out the obvious cases based on scales lying **on** hyperplanes.
  # Now we have to check the hyperplanes that the scales lie above/below.
  svdiff <- signvecX - signvecY

  # If both scales lie off a hyperplane, to be adjacent they should lie in the same direction (diff of 0).
  # But the freer region can step "up" or "down" off a hyperplane of the stricter region (diff of +1 or -1).
  # So what we need to catch are cases where one scale lies above & the other lies below a hyerplane,
  # which shows up as a difference of +2 or -2 in the svdiff.
  difftypes <- unique(abs(svdiff))

  # Want to return -1 if +/-2 is present, else 1 if +/-1 is present, and 0 if identical.
  index <- max(difftypes)
  res <- ((2*index)%%3)-index
  return(res)
}