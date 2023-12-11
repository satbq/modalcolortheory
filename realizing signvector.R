approximate_from_signvector <- function(signvec, ineqmat=NULL, card=NULL, edo=globaledo, rounder=globalrounder) {
  if (is.null(ineqmat)) {
    if (is.null(card)) { warning("Either cardinality or inequality matrix must be specified.") } else {
      ineqmat <- getineqmat(card) }
  }

  res <- qr.solve(ineqmat[, 2:(dim(ineqmat)[2]-1)], signvec) # this assumes T equivalence and a central arrangement
  res <- coord_from_edo(c(0, res))
  return(res)
}

get_relevant_rows <- function(generic_intervals, ineqmat) {
  # Process the generic intervals. They should be indexed such that unisons are 0.
  generic_intervals <- generic_intervals[generic_intervals > 0]
  card <- length(ineqmat[1,]) - 1
  generic_intervals <- abs(sapply(generic_intervals, signed_interval_class, edo=card))

  generics_in_row <- function(row) {
    row <- head(row, -1) # Ignore last column, which doesn't affect generic intervals.
    negative_positions <- which(row < 0)
    postitive_positions <- which(row > 0)
    reference_pitch <- negative_positions[1]
    res <- postitive_positions - reference_pitch
    res <- abs(sapply(res, signed_interval_class, edo=card))
    res <- sort(unique(res))
    return(res)
  }

  row_generics <- apply(ineqmat, 1, generics_in_row)
  check_row <- function(row, generic_interval) generic_interval %in% row
  check_generic <- function(generic_interval) {
    list_of_checks <- lapply(row_generics, check_row, generic_interval=generic_interval)
    return(which(unlist(list_of_checks)==TRUE))
  }

  res <- sapply(generic_intervals, check_generic)
  res <- sort(unique(as.vector(res)))
  return(res)
}

step_signvector <- function(set, edo=globaledo, rounder=globalrounder) {
  card <- length(set)
  ineqmat <- getineqmat(card)

  step_rows <- ineqmat[get_relevant_rows(1, ineqmat=ineqmat),]
  step_sv <- signvector(set, ineqmat=step_rows, edo=edo, rounder=rounder)
  return(step_sv)
}

set_from_signvector <- function(signvec, card, nmax=12, reconvert=FALSE, ineqmat=NULL,
                                edo=globaledo, rounder=globalrounder) {
  if (is.null(ineqmat)) { ineqmat <- makeineqmat(card) }
  set_from_word <- approximate_from_signvector(signvec[get_relevant_rows(1, ineqmat=ineqmat)],
                                               ineqmat=ineqmat[get_relevant_rows(1, ineqmat=ineqmat),])
  implied_word <- asword(set_from_word, edo=edo, rounder=rounder)

  # if (!isTRUE(all.equal(signvec[get_relevant_rows(1, ineqmat=ineqmat)],
  #                       signvector(set_from_word,
  #                                  ineqmat=ineqmat[get_relevant_rows(1, ineqmat=ineqmat)])))) {
  #   warning("Failed to identify implied step word.")
  #   return(rep(NA, card))
  # }

  # From here we reuse code from quantize_color:
  letters <- sort(unique(implied_word), decreasing=FALSE)

  startedo <- sum(implied_word)
  current_set <- cumsum(c(0, implied_word))[1:card]

  if (isTRUE(all.equal(signvector(current_set, ineqmat=ineqmat, edo=startedo, rounder=rounder), signvec))) {
    result_list <- list(set=current_set, edo=startedo)
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
      res <- implied_word

    for (j in seq_along(letters)) {
      res <- replace(res, which(implied_word==letters[j]), newletters[j])
    }

    current_edo <- sum(res)
    current_set <- cumsum(c(0,res))[1:card]

    current_signvec <- signvector(current_set, ineqmat=ineqmat, edo=current_edo, rounder=rounder)

    if (isTRUE(all.equal(current_signvec, signvec))) {
          result_list <- list(set=current_set, edo=current_edo)
          if (reconvert==TRUE) {
            return(convert(result_list$"set", result_list$"edo", edo))
          } else {
            return(result_list)
          }
    }
  }
  return(rep(NA,card))

}