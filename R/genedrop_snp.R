#' Conduct a single Gene-Drop Simulation for a biallelic locus.
#'
#' Conduct a single Gene-Drop Simulation for a biallelic locus.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param genotype vector. Genotypes corresponding to id.
#' @param nsim integer. Number of genedrop simulations to run.
#' @param n_founder_cohorts integer. The number of cohorts at the top of the
#'   pedigree that will sample from the true allele frequencies (these are
#'   defined as "sampled"). All cohorts following these ones are "simulated" and
#'   are used for comparisons of changes in allele frequency.
#' @param fix_founders logical. Default = TRUE. Determines whether individuals
#'   in founder cohorts should be given their true recorded genotypes (if
#'   known). For individuals with no known genotype, their genotypes are sampled
#'   based on the observed cohort allele frequency. If FALSE, then all IDs are
#'   sampled based on the cohort allele frequencies.
#' @param resample_offspring logical. Default = FALSE. If FALSE, the same
#'   pedigree structure as the observed pedigree is used. If TRUE, then
#'   offspring are resampled across parents in each cohort. This is to remove
#'   any potential signal where prolific individuals tend to have prolific
#'   offspring, but will also mean that pedigrees are not directly comparable.
#' @param verbose logical. Default = TRUE. Output the progress of the run.
#' @param interval integer. Default 100. Output progress every 100 simulations.
#' @import dplyr
#' @import kinship2
#' @import magrittr
#' @export
#


genedrop_snp <- function(id,
                         mother,
                         father,
                         cohort = NULL,
                         genotype,
                         nsim,
                         n_founder_cohorts = 1,
                         fix_founders = T,
                         verbose = T,
                         interval = 100,
                         resample_offspring = F){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Check the data and obtain ped object
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ped <- check_data(id, mother, father, cohort, genotype)

  rm(id, mother, father, cohort, genotype)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Get the observed population frequency information
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Summarise the genotype counts per cohort

  x <- table(ped$cohort, ped$genotype, useNA = "always")
  x <- matrix(x, ncol = ncol(x), dimnames = dimnames(x))
  if (any(is.na(row.names(x)))) x <- x[-which(is.na(row.names(x))),]

  x <- cbind(data.frame(cohort = row.names(x)), x)
  names(x)[ncol(x)]<- "NA"

  x$GenoCount <- rowSums(x[,2:(ncol(x)-1)])
  x$FullCount <- rowSums(x[,2:(ncol(x)-1)])
  x$PropGenotyped <- x$GenoCount / x$FullCount

  x$cohort <- as.character(x$cohort)

  #~~ Add any missing genotype columns

  if (is.null(x$`0`)) x$`0` <- 0
  if (is.null(x$`1`)) x$`1` <- 0
  if (is.null(x$`2`)) x$`2` <- 0

  x$p <- (x$`0` + 0.5*(x$`1`)) / x$GenoCount

  # Find cohorts with no representation and throw error

  badcohorts <- x$cohort[which(is.na(x[,2]))]
  if(length(badcohorts) > 0){
    stop(paste0("Cohorts ", paste(badcohorts, collapse = ", "), " have no genotyped individuals."))
  }

  #~~ Create a list to save results

  ped.hold <- ped

  sim.list <- list()

  #~~ Run the simulations

  for(simulation in 1:nsim){

    if (verbose){
      if (simulation %in% seq(1, nsim, interval)){
        message(paste0("Running simulation ", simulation, " of ", nsim, "."))
      }
    }

    if(resample_offspring){
      ped <- resample_offspring_func(ped.hold)
    } else {
      ped <- ped.hold
    }

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # Sample the genotypes in the founder cohorts
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    # #~~ If founders are fixed, pull out the IDs in founder cohorts where genotype is known.
    #
    # if (fix_founders == T){
    #   fixed_founders <- subset(ped,
    #                            cohort %in% x$cohort[1:n_founder_cohorts] &
    #                              !is.na(genotype))$ID
    # } else {
    #   fixed_founders <- NULL
    # }
    #
    #~~ index the pedigree

    row.names(ped) <- ped$ID

    #~~ Create columns for parentally inherited alleles and add alleles for all
    #   individuals that are in the founder cohorts.

    ped$Mum.Allele <- NA
    ped$Dad.Allele <- NA

    if(fix_founders){

      # Generate alleles
      ped$Mum.Allele <- ifelse(ped$genotype %in% 0:1, 0, ifelse(is.na(ped$genotype), NA, 1))
      ped$Dad.Allele <- ifelse(ped$genotype %in% 0  , 0, ifelse(is.na(ped$genotype), NA, 1))

      # Blank out anything that is not in the founder cohorts
      ped$Mum.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA
      ped$Dad.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA

      # Blank out founders that have to inherit an allele from a parent
      ped$Mum.Allele[which(!is.na(ped$MOTHER))] <- NA
      ped$Dad.Allele[which(!is.na(ped$FATHER))] <- NA

    }

    # Convert to character
    ped$ID <- as.character(ped$ID)
    ped$MOTHER <- as.character(ped$MOTHER)
    ped$FATHER <- as.character(ped$FATHER)




    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    #~~ Sample the founders

    for(h in 1:n_founder_cohorts){

      # Mother is known
      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))

      if(length(y1) > 0)  haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      # Father is known
      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))

      if(length(y2) > 0)  haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      # Mother allele is not known - sample from cohort frequency
      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Mum.Allele))

      if(length(y3) > 0){
        haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) ((runif (1) > x$p[h]) + 0L))
      }

      # Father allele is not known - sample from cohort frequency
      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Dad.Allele))

      if(length(y4) > 0){
        haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) ((runif (1) > x$p[h]) + 0L))
      }

      rm(y1, y2, y3, y4)

    }

    #~~ Calculate the cohort frequencies

    cohort.freqs <- data.frame(Sum = tapply(haplo.frame$Mum.Allele, haplo.frame$cohort, sum) +
                                 tapply(haplo.frame$Dad.Allele, haplo.frame$cohort, sum),
                               Count = tapply(haplo.frame$cohort, haplo.frame$cohort, length))

    cohort.freqs$p <- 1 - (cohort.freqs$Sum / (cohort.freqs$Count * 2))

    cohort.freqs$Sum <- NULL
    cohort.freqs$Count <- NULL

    #~~ sample the rest

    for(h in (n_founder_cohorts+1):nrow(x)){

      # Mother is known
      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))

      if(length(y1) > 0)  haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      # Father is known
      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))

      if(length(y2) > 0)  haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      # Estimate the allele frequency
      cohort.freqs$p[h] <- 1 - (sum(haplo.frame$Mum.Allele[y1]) + sum(haplo.frame$Dad.Allele[y2])) / (length(y1) + length(y2))

      if (is.na(cohort.freqs$p[h])){
        stop (paste("Cohort frequency can't be estimated. Problem simulation", simulation, "generation", h))
      }

      # Now sample IDs with missing MOTHER
      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$MOTHER))

      if(length(y3) > 0){
        haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) ((runif (1) > cohort.freqs$p[h]) + 0L))
      }

      # Now sample IDs with missing FATHER
      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$FATHER))

      if(length(y4) > 0){
        haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) ((runif (1) > cohort.freqs$p[h]) + 0L))
      }
      rm(y1, y2, y3, y4)

    }

    head(haplo.frame)

    haplo.frame$MOTHER <- NULL
    haplo.frame$FATHER <- NULL

    haplo.frame$Simulation <- simulation

    sim.list[[simulation]] <- haplo.frame

    rm(cohort.freqs, haplo.frame, h)

  }

  sim.results <- bind_rows(sim.list)

  sim.results$Simulated.Geno <- sim.results$Mum.Allele + sim.results$Dad.Allele

  names(sim.results)[names(sim.results) == "genotype"] <- "True.Geno"

  sim.results

}



