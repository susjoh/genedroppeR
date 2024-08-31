#' Conduct a single Gene-Drop Simulation for a multi-allelic locus.
#'
#' Conduct a single Gene-Drop Simulation for a multi-allelic locus.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year) corresponding to the id.
#' @param genotype vector. Genotypes IDs corresponding to id.
#' @param genotype_delim char. A character denoting the genotype delimited. Default = "".
#' @param nsim integer. Number of genedrop simulations to run.
#' @param n_founder_cohorts integer. The number of cohorts at the top of the
#'   pedigree that will sample from the true allele frequences (these are
#'   defined as "sampled"). All cohorts following these ones are "simulated" and
#'   are used for comparisons of changes in allele frequency.
#' @param fix_founders logical. Default = TRUE. Determines whether individuals in
#'   founder cohorts should be given their true recorded genotypes (if known).
#'   For individuals with no known genotype, their genotypes are sampled based
#'   on the observed cohort allele frequency. If FALSE, then all IDs are sampled
#'   based on the cohort allele frequencies.
#' @param resample_offspring logical. Default = FALSE. If FALSE, the same
#'   pedigree structure as the observed pedigree is used. If TRUE, then
#'   offspring are resampled across parents in each cohort. This is to remove
#'   any potential signal where prolific individuals tend to have prolific
#'   offspring, but will also mean that pedigrees are not directly comparable.
#' @param verbose logical. Output the progress of the run.
#' @param interval int. Default 100. Output progress every 100 simulations.
#' @export


genedrop_multi <- function(id,
                           mother,
                           father,
                           cohort = NULL,
                           genotype,
                           genotype_delim = '',
                           nsim,
                           n_founder_cohorts = 1,
                           fix_founders = TRUE,
                           verbose = TRUE,
                           interval = 100,
                           resample_offspring = F){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Format the data           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ped <- check_data(id, mother, father, cohort, genotype, multiallelic = TRUE)

  rm(id, mother, father, cohort, genotype)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 2. Get the observed population frequency information #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Get the individual allele counts.

  y <- subset(ped, select = c(genotype, cohort))
  y$genotype <- as.character(y$genotype)
  y$Allele1 <- sapply(y$genotype, function(foo)
    strsplit(foo, split = genotype_delim, fixed = T)[[1]][1])
  y$Allele2 <- sapply(y$genotype, function(foo)
    strsplit(foo, split = genotype_delim, fixed = T)[[1]][2])

  y <- melt(y, id.vars = c("genotype", "cohort"))

  x.allele <- sort(unique(y$value))

  x <- table(y$cohort, y$value, useNA = "always")
  x <- matrix(x, ncol = ncol(x), dimnames = dimnames(x))
  if (any(is.na(row.names(x)))) x <- x[-which(is.na(row.names(x))),]


  x <- cbind(data.frame(cohort = row.names(x)), x)
  x$GenoCount <- rowSums(x[,2:(ncol(x)-1)])
  x$FullCount <- rowSums(x[,2:(ncol(x)-1)])

  for(i in x.allele) x[,i] <- x[,i]/x$GenoCount


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 3. Sample the genotypes in the founder cohorts     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ If founders are fixed, determine which IDs are founders!

  if (fix_founders == T){

    fixed_founders <- subset(ped, cohort %in% x$cohort[1:n_founder_cohorts] & !is.na(genotype))$ID

  } else {

    fixed_founders <- NULL

  }

  #~~ index the pedigree

  row.names(ped) <- ped$ID

  #~~ Create columns for parentally inherited alleles and add some for the founders

  ped$Mum.Allele <- NA
  ped$Dad.Allele <- NA

  ped$Mum.Allele <- sapply(ped$genotype, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][1])
  ped$Dad.Allele <- sapply(ped$genotype, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][2])

  ped$Mum.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA
  ped$Dad.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA

  ped$ID <- as.character(ped$ID)
  ped$MOTHER <- as.character(ped$MOTHER)
  ped$FATHER <- as.character(ped$FATHER)


  #~~ Create a list to save results

  ped.hold <- ped


  sim.list <- list()

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

    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    #~~ Sample the founders

    for(h in 1:n_founder_cohorts){

      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER) & !haplo.frame$ID %in% fixed_founders)

      if(length(y1) > 0) haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER) & !haplo.frame$ID %in% fixed_founders)

      if(length(y2) > 0) haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Mum.Allele))

      if(length(y3) > 0) haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) sample(x.allele, size = 1, prob = x[h, x.allele]))

      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Dad.Allele))

      if(length(y4) > 0) haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) sample(x.allele, size = 1, prob = x[h, x.allele]))

      rm(y1, y2, y3, y4)

    }

    #~~ sample the rest

    for(h in (n_founder_cohorts+1):nrow(x)){

      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))

      if(length(y1) > 0) haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))

      if(length(y2) > 0) haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                          1,
                                          function(y) y[((runif (1) > 0.5) + 1L)])

      #~~ Get allele frequencies

      temp.freq <- data.frame(table(c(haplo.frame$Mum.Allele[y1], haplo.frame$Dad.Allele[y2])))
      temp.freq$Freq <- temp.freq$Freq / sum(temp.freq$Freq)
      temp.freq$Var1 <- as.character(temp.freq$Var1)

      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$MOTHER))

      if(length(y3) > 0)  haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$FATHER))

      if(length(y4) > 0) haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      rm(y1, y2, y3, y4)
    }

    haplo.frame$MOTHER <- NULL
    haplo.frame$FATHER <- NULL

    haplo.frame$Simulation <- simulation

    sim.list[[simulation]] <- haplo.frame

    rm(haplo.frame, h)

  }

  sim.results <- do.call(rbind, sim.list)

  sim.results$Simulated.Geno <- paste0(sim.results$Mum.Allele, genotype_delim, sim.results$Dad.Allele)

  names(sim.results)[names(sim.results) == "genotype"] <- "True.Geno"

  sim.results

}


