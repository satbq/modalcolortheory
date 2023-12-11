approximate_from_signvector <- function(signvec, ineqmat=NULL, card=NULL, edo=globaledo, rounder=globalrounder) {
  if (is.null(ineqmat)) { ineqmat <- getineqmat(card) }
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

step_signvector <- function(set, ineqmat=NULL, edo=globaledo, rounder=globalrounder) {
  card <- length(set)
  if (is.null(ineqmat)) { ineqmat <- getineqmat(card) }

  step_rows <- ineqmat[get_relevant_rows(1, ineqmat=ineqmat),]
  step_sv <- signvector(set, ineqmat=step_rows, edo=edo, rounder=rounder)
  return(step_sv)
}

quickcheck <- function(set) {
  card <- length(set)
  newset <- approximate_from_signvector(step_signvector(set), ineqmat=getineqmat(card)[get_relevant_rows(1, getineqmat(card)),])
  return(isTRUE(all.equal(asword(set),asword(newset))))
}