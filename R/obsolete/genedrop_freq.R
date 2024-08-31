#' Conduct a single Gene-Drop Simulation with founder allele frequencies.
#'
#' Conduct a simple gene-drop analysis at a single locus, defining alleles and
#' allele frequencies in the founder cohort.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param genotype_delim char. A character denoting the genotype delimited.
#'   Default = "".
#' @param nsim integer. Number of genedrop simulations to run.
#' @param verbose logical. Output the progress of the run.
#' @param interval int. Default 100. Output progress every 100 simulations.


genedrop_freq <- function(id,
                          mother,
                          father,
                          cohort = NULL,
                          allele_ids,
                          founder_allele_freqs,
                          nsim,
                          genotype_delim = "",
                          verbose = T,
                          interval = 100){


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Format the data           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Check that there are no duplicate IDs.

  if(any(as.numeric(names(table(table(id)))) > 1)) stop("Duplicated values in id")

  #~~ If cohort is provided, then check there are no NA's.

  if(!is.null(cohort)) if(any(is.na(cohort))) message("NAs present in cohort - these individuals will be removed.")


  #~~ Make a ped object.

  ped <- data.frame(ID     = id,
                    MOTHER = mother,
                    FATHER = father)

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA

  #~~ Add cohort information to ped.

  if(!is.null(cohort)) {

    ped$cohort <- cohort

    #~~ remove anything not cohorted

    ped <- subset(ped, !is.na(cohort))

    #~~ Check that parents and offspring are not in the same cohort.

    x <- subset(ped, select = c(ID, MOTHER, FATHER, cohort))
    x1 <- subset(ped, select = c(ID, cohort))
    names(x1) <- c("MOTHER", "Mum.cohort")
    suppressMessages(x <- left_join(x, x1))
    names(x1) <- c("FATHER", "Dad.cohort")
    suppressMessages(x <- left_join(x, x1))

    badids <- c(which(x$cohort == x$Mum.cohort),which(x$cohort == x$Dad.cohort))
    if(length(badids) > 0){
      print(x[badids,])
      stop("Some parents and offspring are in the same cohort: see output for problem lines.")
    }

    rm(x, x1, badids)

  } else {

    if(verbose) message("No cohorts defined: cohorts calculated as individual depth in the pedigree.")
    ped$cohort <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)


  }


  #~~ Get rid of parents that aren't in the IDs.

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 3. Sample the genotypes in the founder cohorts     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  row.names(ped) <- ped$ID

  #~~ Create columns for parentally inherited alleles and add some for the founders

  ped$Mum.Allele <- NA
  ped$Dad.Allele <- NA

  ped$ID <- as.character(ped$ID)
  ped$MOTHER <- as.character(ped$MOTHER)
  ped$FATHER <- as.character(ped$FATHER)


  #~~ Create a list to save results

  sim.list <- list()

  for(simulation in 1:nsim){

    if(verbose) if(simulation %in% seq(1, nsim, interval)) message(paste0("Running simulation ", simulation, " of ", nsim, "."))
    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    x <- summary_cohort(id = ped$ID, mother = ped$MOTHER, father = ped$FATHER, cohort = ped$cohort)

    #~~ Sample the founders

    y1 <- which(haplo.frame$cohort %in% x$cohort[1])

    haplo.frame$Mum.Allele[y1] <- sapply(y1, function(y) sample(allele_ids, size = 1, prob = founder_allele_freqs))
    haplo.frame$Dad.Allele[y1] <- sapply(y1, function(y) sample(allele_ids, size = 1, prob = founder_allele_freqs))

    rm(y1)

    x <- subset(x, !is.na(cohort))

    for(h in 2:nrow(x)){

      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))
      haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")], 1, function(y) y[((runif(1) > 0.5) + 1L)])

      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))
      haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")], 1, function(y) y[((runif(1) > 0.5) + 1L)])

      #~~ Get allele frequencies

      temp.freq <- data.frame(table(c(haplo.frame$Mum.Allele[y1], haplo.frame$Dad.Allele[y2])))
      temp.freq$Freq <- temp.freq$Freq/sum(temp.freq$Freq)
      temp.freq$Var1 <- as.character(temp.freq$Var1)

      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$MOTHER))
      haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$FATHER))
      haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      rm(y1, y2, y3, y4)
    }

    head(haplo.frame)

    haplo.frame$MOTHER <- NULL
    haplo.frame$FATHER <- NULL

    haplo.frame$Simulation <- simulation

    sim.list[[simulation]] <- haplo.frame

    rm(haplo.frame, h)

  }

  sim.results <- do.call(rbind, sim.list)
  sim.results$Simulated.Geno <- paste0(sim.results$Mum.Allele, genotype_delim, sim.results$Dad.Allele)

  sim.results

}

