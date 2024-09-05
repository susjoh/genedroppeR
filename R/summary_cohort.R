#' Summarise cohorts.
#'
#' Summarise cohorts by number of founders, proportion of individuals genotyped
#' and/or allele frequency information at a given locus.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param genotype vector. Genotypes IDs corresponding to id.
#' @param genotype_delim char. A character denoting the genotype delimited.
#'   Default = ''.
#' @examples
#'
#' data(unicorn)
#' unicorn_summary <- summary_cohort(
#'   id = unicorn$id,
#'   mother = unicorn$mother,
#'   father = unicorn$father,
#'   cohort = unicorn$cohort,
#'   genotype = unicorn$Horns
#' )
#' plot_genedrop_cohort(unicorn_summary)
#'
#' @export


summary_cohort <- function(id,
                           mother = NULL,
                           father = NULL,
                           cohort = NULL,
                           genotype = NULL,
                           genotype_delim = "") {
  # ~~ Check that there are no duplicate IDs.

  if (any(as.numeric(names(table(table(id)))) > 1)) {
    stop("Duplicated values in id")
  }

  # ~~ If genotype is provided, then check the locus is not monomorphic.

  if (!is.null(genotype) & length(table(genotype)) == 1) {
    stop("Locus is monomorphic")
  }


  # ~~ If genotype is numeric, then only accept if 0, 1, 2

  if (!is.null(genotype) & is.numeric(genotype)) {
    if (any(!na.omit(genotype) %in% 0:2)) {
      stop("Dosage should be coded as 0, 1, 2.")
    }

    genotype[genotype == 0] <- "AA"
    genotype[genotype == 1] <- "AB"
    genotype[genotype == 2] <- "BB"

    message("Locus summary AA, AB, BB corresponds to 0, 1, 2")
  }

  # ~~ If no parents are defined, output a message...

  if (is.null(mother) & is.null(father)) {
    message("No parents defined: proportion of founders cannot be calculated.")
  }

  # ~~ Make a ped object.

  if (is.null(mother)) mother <- rep(NA, times = length(id))
  if (is.null(father)) father <- rep(NA, times = length(id))


  ped <- data.frame(
    ID = id,
    MOTHER = mother,
    FATHER = father
  )

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA


  # ~~ Add cohort information to ped.

  if (!is.null(cohort)) {
    ped$cohort <- cohort
  } else {
    ped$cohort <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)
  }


  # ~~ Add the genotype information.


  if (is.null(genotype)) {
    ped2 <- table(ped$cohort, useNA = "always")
    ped2 <- matrix(ped2, ncol = 1, dimnames = list(names(ped2), "FullCount"))
    if (any(is.na(row.names(ped2)))) row.names(ped2)[which(is.na(row.names(ped2)))] <- "missing"
    ped2 <- data.frame(ped2)
    ped2$cohort <- row.names(ped2)
  } else {
    ped$genotype <- as.character(genotype)
    ped$Allele1 <- sapply(ped$genotype, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][1])
    ped$Allele2 <- sapply(ped$genotype, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][2])

    ped2 <- melt(ped[, c("cohort", "Allele1", "Allele2")], id.vars = "cohort")

    x.allele <- sort(unique(ped2$value))

    ped2 <- table(ped2$cohort, ped2$value, useNA = "always")
    ped2 <- matrix(ped2, ncol = ncol(ped2), dimnames = dimnames(ped2))
    if (any(is.na(row.names(ped2)))) row.names(ped2)[which(is.na(row.names(ped2)))] <- "missing"
    ped2 <- data.frame(ped2)
    ped2 <- cbind(cohort = row.names(ped2), ped2)

    ped2$GenoCount <- rowSums(ped2[, 2:(ncol(ped2) - 1)]) / 2

    ped2$FullCount <- rowSums(ped2[, 2:(ncol(ped2) - 1)]) / 2
    ped2$PropGenotyped <- ped2$GenoCount / ped2$FullCount

    for (i in x.allele) ped2[, i] <- (0.5 * ped2[, i]) / ped2$GenoCount

    ped2$NA. <- NULL
  }


  # ~~ How many founders?

  ped$founder <- kindepth(ped[, 1], ped[, 2], ped[, 3])

  suppressMessages(x1 <- table(ped$cohort, ped$founder == 0) %>% data.frame() %>% dcast(Var1 ~ Var2))

  if (ncol(x1) == 2) {
    x1$V3 <- NA
    x1 <- x1[, c(1, 3, 2)]
  }

  names(x1) <- c("cohort", "NonFounders", "Founders")

  x1$PropFounders <- x1$Founders / (x1$NonFounders + x1$Founders)

  suppressMessages(ped2 <- left_join(ped2, x1))

  # ped2$cohort[ped2$cohort == "missing"] <- NA
  # ped2$cohort <- as.numeric(as.character(ped2$cohort))

  ped2
}
