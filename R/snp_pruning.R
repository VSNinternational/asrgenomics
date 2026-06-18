#' Reduces the number of redundant markers on a molecular matrix M by pruning
#'
#' For a given molecular dataset \eqn{\boldsymbol{M}} (in the format 0, 1 and 2)
#' it produces a reduced molecular matrix by eliminating "redundant"
#' markers using pruning techniques. This function finds and drops some of the
#' SNPs in high linkage disequilibrium (LD).
#'
#' Pruning is recommended as redundancies can affect
#' the quality of matrices used for downstream analyses.
#' The algorithm used is based on the Pearson's correlation between markers
#' as a \emph{proxy} for LD. In the event of a pairwise correlation higher
#' than the selected threshold markers will be eliminated as specified by: call rate,
#' minor allele frequency. In case of tie, one marker will be dropped at random.
#'
#' @param M A matrix with marker data of full form (\eqn{n \times p}), with \eqn{n} individuals
#' and \eqn{p} markers. Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows.
#' If \code{NULL} a dummy map is generated considering a single chromosome and sequential positions
#' for markers. A \code{map} is mandatory if \code{by.chrom = TRUE}, where also option \code{chrom}
#' must also be non-null.
#' @param marker A character indicating the name of the column in data frame \code{map}
#' with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param chrom A character indicating the name of the column in data frame \code{map} with the identification
#' of chromosomes. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param pos A character indicating the name of the column in data frame \code{map} with the identification
#' of marker positions (default = \code{NULL}).
#' @param method A character indicating the method (or algorithm) to be used as reference for
#' identifying redundant markers.
#' The only method currently available is based on correlations (default = \code{"correlation"}).
#' @param criteria A character indicating the criteria to choose which marker to drop
#' from a detected redundant pair.
#' Options are: \code{"callrate"} (the marker with fewer missing values will be kept) and
#' \code{"maf"} (the marker with higher minor allele frequency will be kept) (default = \code{"callrate"}).
#' @param pruning.thr A threshold value to identify redundant markers with Pearson's correlation larger than the
#' value provided (default = \code{0.95}).
#' @param by.chrom If TRUE the pruning is performed independently by chromosome (default = \code{FALSE}).
#' @param window.n A numeric value with number of markers to consider in each
#' window to perform pruning (default = \code{50}).
#' @param overlap.n A numeric value with number of markers to overlap between consecutive windows
#' (default = \code{5}).
#' @param iterations An integer indicating the number of sequential times the pruning procedure
#' should be executed on remaining markers.
#' If no markers are dropped in a given iteration/run, the algorithm will stop (default = \code{10}).
#' @param seed An integer to be used as seed for reproducibility. In case the criteria has the
#' same values for a given pair of markers, one will be dropped at random (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details Filtering markers (\link{qc.filtering}) is of high relevance before pruning.
#' Poor quality markers (\emph{e.g.}, monomorphic markers) may prevent correlations from being
#' calculated and may affect eliminations.
#'
#' @return
#' \itemize{
#'  \item{\code{Mpruned}: a matrix containing the pruned marker \emph{M} matrix.}
#'  \item{\code{map}: an data frame containing the pruned map.}
#' }
#'
#' @export
#'
#' @examples
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.20, ind.callrate = 0.20,
#'  Fis = 1, heterozygosity = 0.98,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Prune correlations > 0.9.
#' Mpr <- snp.pruning(
#'  M = M.clean, pruning.thr = 0.90,
#'  by.chrom = FALSE, window.n = 40, overlap.n = 10)
#' head(Mpr$map)
#' Mpr$Mpruned[1:5, 1:5]
#'
snp.pruning <- function(
  M = NULL,
  map = NULL,
  marker = NULL,
  chrom = NULL,
  pos = NULL,
  method = c('correlation'),
  criteria = c("callrate", "maf"),
  pruning.thr = 0.95,
  by.chrom = FALSE,
  window.n = 50,
  overlap.n = 5,
  iterations = 10,
  seed = NULL,
  message = TRUE
) {
  check.data_(data_ = "M", class_ = "matrix")
  maf <- maf(M = M)
  callrate <- callrate(M = M, margin = "col")
  if (any(callrate == 0)) {
    stop(
      "There are markers will all samples missing. Please use qc.filtering() before pruning."
    )
  }
  if (any(maf == 0)) {
    stop(
      "There are markers with minor allele frequency equal to 0. Please use qc.filtering() before pruning."
    )
  }
  if (is.null(map)) {
    map <- dummy.map_(colnames(M))
    marker <- "marker"
    chrom <- "chrom"
    pos <- "pos"
  } else {
    check.data_(data_ = "map", class_ = "data.frame")
    if (is.null(marker)) {
      stop("The \'marker' option must be specified if \'map' is provided.")
    }
    if (is.null(chrom)) {
      stop("The \'chrom' option must be specified if \'map' is provided.")
    }
    map.name.hit <- c(marker, chrom, pos) %in% names(map)
    if (!all(map.name.hit)) {
      stop(
        "Value provided to argument \'",
        c("marker", "chrom", "pos")[!map.name.hit],
        "' does not correspond to a variable in
           data frame \'map'."
      )
    }
    if (!identical(as.character(map[[marker]]), colnames(M))) {
      stop(
        "map[[marker]] and colnames(M) must be identical. Please check input."
      )
    }
  }
  method <- match.arg(method)
  criteria <- match.arg(criteria)
  if (pruning.thr <= 0 | pruning.thr > 1) {
    stop("The condition for pruning.thr is between 0 and 1.")
  }
  check.logical_(arg_ = "by.chrom")
  if (window.n <= 1) {
    stop("The \'window.n' argument should have an integer larger than 1.")
  }
  if (overlap.n <= 0) {
    stop(
      "The \'overlap.n' argument should have an integer larger or eqaul than 0."
    )
  }
  if (overlap.n >= window.n) {
    stop(
      "The \'overlap.n' argument should be lower than the \'window.n' argument."
    )
  }
  if (iterations <= 0) {
    stop("The \'iterations' argument should have an positive integer.")
  }
  check.logical_(arg_ = "message")
  if (!is.null(seed)) {
    set.seed(seed = seed)
  }
  if (criteria == "maf") {
    map$criteria <- maf
  }
  if (criteria == "callrate") {
    map$criteria <- callrate
  }
  map$sel <- 1
  map$tiebreak <- sample.int(nrow(map))
  rm(maf, callrate)
  marker.drop <- function(curr.set.index = NULL) {
    init.set.pos <- sets[curr.set.index]
    if (n.sets == 0) {
      window.M <- cur.M
    } else if (curr.set.index == n.sets) {
      window.M <- cur.M[, sets[curr.set.index]:ncol(cur.M)]
    } else {
      window.M <- cur.M[,
        sets[curr.set.index]:(sets[curr.set.index + 1] + overlap.n - 1)
      ]
    }
    C.sparse <- suppressWarnings(cor(window.M, use = 'pairwise.complete.obs'))
    C.sparse[is.na(C.sparse)] <- 0
    C.sparse <- as.data.table(full2sparse(C.sparse))
    C.sparse[, Value := abs(Value)]
    C.sparse <- C.sparse[Row != Col, ]
    setorder(C.sparse, -Value)
    rm.pos <- c()
    while (
      nrow(C.sparse) > 0L &&
        C.sparse$Value[1] >= pruning.thr
    ) {
      row.pos <- C.sparse[1, Row] + init.set.pos - 1
      col.pos <- C.sparse[1, Col] + init.set.pos - 1
      if (cur.map$criteria[row.pos] == cur.map$criteria[col.pos]) {
        if (cur.map$tiebreak[row.pos] > cur.map$tiebreak[col.pos]) {
          C.sparse <- C.sparse[Col != Col[1] & Row != Col[1], ]
          rm.pos <- append(rm.pos, col.pos)
        } else {
          C.sparse <- C.sparse[Col != Row[1] & Row != Row[1], ]
          rm.pos <- append(rm.pos, row.pos)
        }
      } else if (cur.map$criteria[row.pos] > cur.map$criteria[col.pos]) {
        C.sparse <- C.sparse[Col != Col[1] & Row != Col[1], ]
        rm.pos <- append(rm.pos, col.pos)
      } else {
        C.sparse <- C.sparse[Col != Row[1] & Row != Row[1], ]
        rm.pos <- append(rm.pos, row.pos)
      }
    }
    return(rm.pos)
  }
  if (message) {
    original.n.markers <- ncol(M)
    if (by.chrom) {
      original.n.markers.chrom <- table(map$chrom)
    }
  }
  if (message) {
    message(blue("\nInitiating pruning procedure."))
    message(
      "Initial marker matrix M contains ",
      nrow(M),
      " individuals and ",
      ncol(M),
      " markers."
    )
  }
  if (by.chrom) {
    chrom.range <- unique(map[[chrom]])
    if (message) {
      message("Requesting pruning by chromosome.")
    }
  } else {
    chrom.range <- 1
    if (message) {
      message("Requesting pruning without chromosome indexing.")
    }
  }
  iter.range <- 1:iterations
  for (cur.chrom in chrom.range) {
    if (length(chrom.range) > 1 & message) {
      message(paste0("Chromosome: ", cur.chrom))
    }
    if (by.chrom) {
      split.index <- map[, chrom] == cur.chrom
      cur.map <- map[split.index, ]
      cur.M <- M[, split.index, drop = FALSE]
      map <- map[!split.index, ]
      M <- M[, !split.index, drop = FALSE]
    } else {
      cur.map <- map
      map <- NULL
      cur.M <- M
      M <- NULL
    }
    for (iter in iter.range) {
      if (message) {
        message("  Iteration: ", iter)
      }
      step <- window.n - overlap.n
      sets <- seq(1, ncol(cur.M), step)
      n.sets <- length(sets) - 1
      sets.range <- if (n.sets == 0L) 1L else seq_len(n.sets)
      drop.pos <- unique(unlist(lapply(X = sets.range, FUN = marker.drop)))
      if (length(drop.pos) == 0L) {
        break
      }
      if (length(drop.pos) >= ncol(cur.M)) {
        stop("Internal pruning error: all markers were selected for removal.")
      }
      cur.map$sel[drop.pos] <- 0
      cur.M <- cur.M[, cur.map$sel == 1, drop = FALSE]
      cur.map <- cur.map[cur.map$sel == 1, , drop = FALSE]
    }
    map <- rbind(map, cur.map)
    M <- cbind(M, cur.M)
  }
  if (message) {
    message(
      "\nFinal pruned marker matrix M contains ",
      nrow(M),
      " individuals and ",
      ncol(M),
      " markers."
    )
    message(
      "A total of ",
      original.n.markers - ncol(M),
      " markers were pruned."
    )
    if (by.chrom) {
      message(paste0(
        "A total of ",
        table(map$chrom),
        " markers were kept in chromosome ",
        names(table(map$chrom)),
        ".",
        collapse = "\n"
      ))
      message(paste0(
        "A total of ",
        original.n.markers.chrom - table(map$chrom),
        " markers were pruned from chromosome ",
        names(table(map$chrom)),
        ".",
        collapse = "\n"
      ))
    }
    message(
      "Range of minor allele frequency after pruning: ",
      paste0(round(range(maf(M = M)), 2), collapse = " ~ ")
    )
    message(
      "Range of marker call rate after pruning: ",
      paste0(round(range(callrate(M = M, margin = "col")), 2), collapse = " ~ ")
    )
    message(
      "Range of individual call rate after pruning: ",
      paste0(round(range(callrate(M = M, margin = "row")), 2), collapse = " ~ ")
    )
  }
  map <- map[, !names(map) %in% c("criteria", "sel", "tiebreak")]
  return(list(map = map, Mpruned = M))
}
