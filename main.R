globaledo <- 12
tiny <- 1e-10
globalrounder <- -log(tiny)/log(10)
fortenums <- readRDS("fortenums.rds")
ineqmats <- readRDS("ineqmats.rds")
representative_signvectors <- readRDS("representative_signvectors.rds")
representative_scales <- readRDS("representative_scales.rds")

# Useful Scalar Constants. Note that these don't automatically update if you redefine globaledo.
limma = globaledo * log(256/243)/log(2)
just_st = globaledo * log(16/15)/log(2)
apotome = globaledo * log(2187/2048)/log(2)
just_wt = globaledo * log(9/8)/log(2)
just_min3 = globaledo * log(6/5)/log(2)
just_maj3 = globaledo * log(5/4)/log(2)
just_p4 = globaledo * log(4/3)/log(2)
just_p5 = globaledo * log(3/2)/log(2)

pyth_comma = (12 * just_fifth) %% globaledo
pyth_maj3 = (4 * just_p5) %% globaledo
syntonic_comma = pyth_maj3 - just_maj3
meantone_fifth = function(frac=1/4) just_fifth - (syntonic_comma * frac)


# Main Functions

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
  res <- apply(res, 2, startzero, edo, sorted=FALSE)
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
  upset <- startzero(tnprime(set, edo), edo)
  downset <- startzero(tnprime(tni(set, 0, edo), edo), edo)
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

# "maximally even MOS"
makeMEMOS <- function(card, edo=globaledo, floor=TRUE) {
  if (floor==TRUE) {
    res <- primeform(floor(edoo(card, edo)), edo)
  } else {
    res <- primeform(round(edoo(card, edo), digits=0), edo)
  }

  return(res)
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

makesubmat <- function(card) {
  submat <- t(outer(0:(card-1), 0:(card-1), paste, sep="-"))

  rotaterow <- function(n,mat) rotate(mat[n,], n-1)
  submat <- t(sapply(1:card, rotaterow,mat=submat))

  submat <- submat[,card:1]
  submat[lower.tri(submat)] <- mapply(paste0, "w+", submat[lower.tri(submat)])
  submat <- submat[,card:1]

  return(submat)
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

realize_setword <- function(setword, edo=globaledo) {
  set <- cumsum(setword)
  wordedo <- set[length(set)]
  set <- convert(sort(set %% wordedo), wordedo, edo)
  return(set)
}

iswellformed <- function(set, setword=NULL, allowdegen=FALSE, edo=globaledo, rounder=globalrounder) {
  if ( is.null(set) ) { set <- realize_setword(setword, edo) }
  speccount <- spectrumcount(set, edo, rounder)[-1]
  uniques <- unique(speccount)
  if (toString(uniques)=="2") { return(TRUE) }
  if (toString(uniques)=="1") { return(as.logical(allowdegen)) }
  if (toString(sort(uniques))=="1, 2") { return(as.logical(allowdegen)) }

  return(FALSE)
}

equivocate <- function(setword, lowerbound, windowsize) {
  highest <- max(setword)
  toMatch <- lowerbound:(lowerbound+(windowsize-1))
  toMatch <- unique(((toMatch-1)%%highest)+1)
  replacement_positions <- which(setword %in% toMatch)
  result <- replace(setword, replacement_positions, 1)
  result <- replace(result, -replacement_positions, 2)
  return(result)
}

isgwf <- function(set, setword=NULL,allowdegen=FALSE,edo=globaledo,rounder=globalrounder) {
# Note that this requires that the "letters" of a setword are consecutive integers,
# such that max(letters) == len(unique(letters)). That is, a word on 3 letters should have the letters 1, 2, and 3.
  if ( is.null(setword) ) { setword <- asword(set, edo, rounder) }
  if (anyNA(setword)) { return(FALSE) }

  highest <- max(setword)
  equiv_parameters <- expand.grid(1:highest, 1:(highest-1))

  equiv_wrap <- function(params, setword) equivocate(setword, params[1], params[2])
  reduced_words <- apply(equiv_parameters, 1, equiv_wrap, setword=setword)

  iswf_wrap <- function(setword, allowdegen, edo, rounder)  {
    return(iswellformed(NULL, setword, allowdegen, edo, rounder))
  }

  tests <- apply(reduced_words,2, iswf_wrap, allowdegen=allowdegen, edo=edo, rounder=rounder)

  return(as.logical(prod(tests)))
}

OPTC_test <- function(set, edo=globaledo, rounder=globalrounder, single_answer=TRUE) {
  basically_zero <- 10^(-rounder)
  step_sizes <- sim(set)[2,]

  satisfies_O <- max(set) < edo
  satisfies_P <- isTRUE( all.equal(set, sort(set), tolerance=basically_zero) )
  satisfies_T <- abs(set[1]) < basically_zero
  satisfies_C <- min(step_sizes) > basically_zero

  if (single_answer == FALSE) { return(c(satisfies_O, satisfies_P, satisfies_T, satisfies_C)) }
  return(satisfies_O && satisfies_P && satisfies_T && satisfies_C)
}

quantize_color <- function(set,nmax=12,reconvert=FALSE,edo=globaledo,rounder=globalrounder) {
  card <- length(set)
  signvec <- signvector(set,edo=edo,rounder=rounder)

  word <- asword(set, edo, rounder)
  letters <- sort(unique(word),decreasing=FALSE)

  startedo <- sum(word)

  current_set <- cumsum(c(0,word))[1:card]

  if (isTRUE(all.equal(signvector(current_set, edo, rounder), signvec))) {
    result_list <- list(set=curset, edo=startedo)
    if (reconvert==TRUE) {
      return(convert(result_list$"set", result_list$"edo", edo))
    } else {
      return(result_list)
    }
  }

  options <- combn(nmax,length(letters))
  stop <- dim(options)[2]

  for (i in 1:stop) {
      newletters <- options[,i]
      res <- word

    for (j in seq_along(letters)) {
      res <- replace(res, which(word==letters[j]), newletters[j])
    }

    current_edo <- sum(res)
    current_set <- cumsum(c(0,res))[1:card]

    current_signvec <- signvector(current_set, edo=current_edo, rounder=rounder)

    if (isTRUE(all.equal(current_signvec, signvec))) {
          result_list <- list(set=current_set, edo=current_edo)
          if (reconvert==TRUE) {
            return(convert(result_list$"set", result_list$"edo", edo))
          } else {
            return(result_list)
          }
    }
  }
  return(NA)
}

colornum <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder,
                     signvector_list=representative_signvectors) {
  if (evenness(set, edo) < 10^(-rounder) ) { return(0) }

  card <- length(set)
  signvec <- toString(signvector(set, ineqmat, edo, rounder))

  return(which(signvector_list[[card]]==signvec))
}


# Scala-Related Functions

writeSCL <- function(x, filename, period=2, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  # Period defined as a frequency ratio (i.e. 2 for octave-repeading scales)
  periodCents <- 1200 * log(period)/log(2)

  if (missing(filename)) {
    filename <- deparse(substitute(x))
    filename <- paste0(filename, ".scl")
  }

  ###Comment with the filename

  line0 <- paste0("! ",filename)

  ###The scale description
  card <- length(x)
  freedom <- howfree(x)
  svzeroes <- countsvzeroes(x, ineqmat=ineqmat, edo=edo, rounder=rounder)
  scaleratio <- round(ratio(x, edo=edo, rounder=rounder),3)
  how_even <- round(evenness(x, edo=edo), 3)

  nameColor <- FALSE
  if (exists(deparse(substitute(representative_signvectors)))) {
    color <- colornum(x,ineqmat=ineqmat, edo=edo, rounder=rounder,
                     signvector_list=representative_signvectors)
    nameColor <- TRUE
  }


  if (freedom==1) { degree <- "degree" } else { degree <- "degrees" }
  if (svzeroes==1) { zeroes <- "hyperplane" } else { zeroes <- "hyperplanes" }

  if (nameColor==TRUE) {
    line1 <- paste0("A ", card, "-note scale of color ", color, " with ", freedom, " ", degree, " of freedom. Period: ", round(periodCents,3), " cents. Ratio: ", scaleratio, "; evenness: ", how_even, " (wrt. ", edo, " steps to period); on ", svzeroes, " ", zeroes, ".")
  } else {
    line1 <- paste0("A ", card, "-note scale with ", freedom, " ", degree, " of freedom. Period: ", round(periodCents,3), " cents. Ratio: ", scaleratio, "; evenness: ", how_even, " (wrt. ", edo, " steps to period); on ", svzeroes, " ", zeroes, ".")
  }

  if (how_even < 10^-rounder) { line1 <- paste0("A ", card, "-note equal division of ", round(periodCents, 3), " cents.") }

  ###The data for Scala

  line2 <- paste(card)
  line3 <- "! "
  x <- c(x,edo)
  scaledef <- convert(x,edo,periodCents)[-1]
  scaledef <- format(scaledef,nsmall=1,digits=(rounder+2))

  fileConn <- file(filename)
  writeLines(c(line0,line1,line2,line3,scaledef),fileConn)
  close(fileConn)
}


readSCL <- function(filename, scaleonly=TRUE, edo=globaledo) {
#Note that if scala files aren't defined with enough precision, they may appear to lack structure.
#Check how close they are to hyperplanes to adjust rounder and tiny appropriately.
  contents <- scan(filename,what="character",sep="\n",quiet=TRUE,blank.lines.skip=FALSE)

  removeAfterChar <- function(string,charToRemove) {
    charindex <- as.vector(regexpr(charToRemove,string,fixed=TRUE))
    charindex[charindex ==-1] <- nchar(string[charindex==-1]) + 1
    charindex <- charindex - 1
    res <- do.call(substr,list(x=string,start=rep(0,length(string)),stop=charindex))
    res <- res[res!=""]
    return(res)
  }

  #Fill in description line if blank
  bangindex <- as.vector(grepl("!",contents,fixed=TRUE))
  firstrealline <- contents[which(bangindex==FALSE)[1]]
  if (firstrealline == "") { contents[which(bangindex==FALSE)[1]] <- "blank description" }

  #Remove comments marked by !
  contents <- removeAfterChar(contents,"!")
  contents <- trimws(contents)

  #Read scale length & remove description header
  card <- strtoi(contents[2])
  if (class(card) != "integer") { warning(".scl file not formatted as expected. Scale length not found.") }
  contents <- contents[-(1:2)]

  #Locate the degrees defined as ratios vs. those defined by cents
  ratioindex <- grepl("/",contents,fixed=TRUE)
  centsindex <- grepl(".",contents,fixed=TRUE)

  #Remove whitespace-prefixed comments
  contents <- trimws(contents)
  contents <- removeAfterChar(contents," ")

  #Check for integers (e.g. an octave defined as "2" rather than "2/1")
  integerindex <- !is.na(strtoi(contents))

  if (sum(ratioindex,centsindex,integerindex) != card) { warning(".scl file not formatted as expected. Cents or ratios incorrectly identified.") }

  #Convert the ratios to scale degrees
  ratioToLog <- function(ratioString, edo=edo) {
    dividends <- unlist(strsplit(ratioString, split="/"))
    if ( length(dividends) != 2 ) { warning(".scl file not formatted as expected. Apparent ratio with more or less than 2 arguments.") }
    dividends <- as.numeric(dividends)
    res <- dividends[1]/dividends[2]
    res <- edo * log(res)/log(2)
    return(res)
  }

  #Convert integers to scale degrees
  intToLog <- function(intString, edo=edo) {
    val <- strtoi(intString)
    return(edo * log(val)/log(2))
  }

  #Apply above functions and convert cents to global EDO
  if (sum(ratioindex) > 0) { contents[ratioindex] <- sapply(contents[ratioindex],ratioToLog, edo=edo) }
  if (sum(integerindex) > 0) {contents[integerindex] <- sapply(contents[integerindex],intToLog) }
  if (sum(centsindex) > 0) { contents[centsindex] <- sapply(contents[centsindex],as.numeric) * (edo/1200) }
  contents <- as.numeric(contents)

  #Format & output results
  names(contents) <- NULL
  octave <- contents[length(contents)]
  contents <- c(0,contents[-length(contents)])
  if (scaleonly==TRUE) { return(contents) }
  return(list(scale=contents,length=card,period=octave))
}