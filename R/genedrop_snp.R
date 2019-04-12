#' Conduct a single Gene-Drop Simulation for a biallelic locus
#'
#' Conduct a single Gene-Drop Simulation for a biallelic locus
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year) corresponding to the id.
#' @param genotype vector. Genotypes IDs corresponding to id.
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
#' @param verbose logical. Output the progress of the run.
#' @param interval int. Default 100. Output progress every 100 simulations.
#' @import plyr
#' @import kinship2
#' @import magrittr
#' @export



genedrop_snp <- function(id,
                         mother,
                         father,
                         cohort = NULL,
                         genotype,
                         nsim,
                         n_founder_cohorts = 1,
                         fix_founders = T,
                         verbose = T,
                         interval = 100){

  require(kinship2)
  require(magrittr)
  require(plyr)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Format the data           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Check that there are no duplicate IDs.

  if(any(as.numeric(names(table(table(id)))) > 1)) stop("Duplicated values in id")

  #~~ If cohort is provided, then check there are no NA's.

  if(!is.null(cohort) & any(is.na(cohort))) message("NAs present in cohort - these individuals will be removed.")

  #~~ If genotype is provided, then check the locus is not monomorphic.

  if(!is.null(genotype) & length(table(genotype)) == 1) stop("Locus is monomorphic")

  #~~ If genotype is numeric, then only accept if 0, 1, 2

  if(!is.null(genotype) & is.numeric(genotype)){

    if(any(!na.omit(genotype) %in% 0:2)){
      stop("Dosage is only tolerated for biallelic markers and should be coded as 0, 1, 2.")
    }

  } else {

    templev <- as.factor(genotype) %>% levels %>% sort

    #~~ determine which value is heterozygote if templev is shorter than 3

    if(length(templev) == 2){

      templev2 <- sapply(templev, function(x) length(unique(strsplit(x, split = "")[[1]])))
      if(templev2[1] == max(templev2)) templev <- c("AA", names(templev2))
      if(templev2[1] == templev2[2]) templev <- c(names(templev2)[1], "AB", names(templev2)[2])

      rm(templev2)
    }

    genotype <- as.numeric(factor(genotype, levels = templev)) -1
    rm(templev)
  }

  #~~ Make a ped object.

  ped <- data.frame(ID     = id,
                    MOTHER = mother,
                    FATHER = father)

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA

  #~~ Add the genotype information.

  if(!is.null(genotype)) ped$genotype <- genotype

  #~~ Add cohort information to ped.

  if(!is.null(cohort)) {

    ped$cohort <- cohort

    #~~ remove anything not cohorted

    ped <- subset(ped, !is.na(cohort))

    #~~ Check that parents and offspring are not in the same cohort.

    x <- subset(ped, select = c(ID, MOTHER, FATHER, cohort))
    x1 <- subset(ped, select = c(ID, cohort))
    names(x1) <- c("MOTHER", "Mum.cohort")
    suppressMessages(x <- join(x, x1))
    names(x1) <- c("FATHER", "Dad.cohort")
    suppressMessages(x <- join(x, x1))

    badids <- c(which(x$cohort == x$Mum.cohort),which(x$cohort == x$Dad.cohort))
    if(length(badids) > 0){
      print(x[badids,])
      stop("Some parents and offspring are in the same cohort: see output for problem lines.")
    }

    rm(x, x1, badids)

  } else {

    ped$cohort <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)


  }


  #~~ Get rid of parents that aren't in the IDs.

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 2. Get the observed population frequency information #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Summarise the genotype counts

  x <- table(ped$cohort, ped$genotype, useNA = "always")
  x <- matrix(x, ncol = ncol(x), dimnames = dimnames(x))
  if(any(is.na(row.names(x)))) x <- x[-which(is.na(row.names(x))),]


  x <- cbind(data.frame(cohort = row.names(x)), x)
  names(x)[ncol(x)]<- "NA"

  x$GenoCount <- rowSums(x[,2:(ncol(x)-1)])
  x$FullCount <- rowSums(x[,2:(ncol(x)-1)])
  x$PropGenotyped <- x$GenoCount/x$FullCount

  x$cohort <- as.character(x$cohort)

  x

  #~~ Add any missing genotypes

  if(is.null(x$`0`)) x$`0` <- 0
  if(is.null(x$`1`)) x$`1` <- 0
  if(is.null(x$`2`)) x$`2` <- 0

  x$p <- (x$`0` + 0.5*(x$`1`))/x$GenoCount

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 3. Sample the genotypes in the founder cohorts     #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ If founders are fixed, determine which IDs are founders!

  if(fix_founders == T){
    fixed_founders <- subset(ped, cohort %in% x$cohort[1:n_founder_cohorts] & !is.na(genotype))$ID
  } else {
    fixed_founders <- NULL
  }

  #~~ index the pedigree

  row.names(ped) <- ped$ID

  #~~ Create columns for parentally inherited alleles and add some for the founders

  ped$Mum.Allele <- NA
  ped$Dad.Allele <- NA

  ped$Mum.Allele <- ifelse(ped$genotype %in% 0:1, 0, ifelse(is.na(ped$genotype), NA, 1))
  ped$Dad.Allele <- ifelse(ped$genotype %in% 0  , 0, ifelse(is.na(ped$genotype), NA, 1))

  ped$Mum.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA
  ped$Dad.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA

  ped$ID <- as.character(ped$ID)
  ped$MOTHER <- as.character(ped$MOTHER)
  ped$FATHER <- as.character(ped$FATHER)


  #~~ Create a list to save results

  sim.list <- list()

  for(simulation in 1:nsim){

    if(verbose) if(simulation %in% seq(1, nsim, interval)) message(paste0("Running simulation ", simulation, " of ", nsim, "."))
    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    #~~ Sample the founders

    for(h in 1:n_founder_cohorts){

      y1 <- which(haplo.frame$cohort %in% x$cohort[h] & !haplo.frame$ID %in% fixed_founders)

      haplo.frame$Mum.Allele[y1] <- sapply(y1, function(y) ((runif(1) > x$p[h]) + 0L))
      haplo.frame$Dad.Allele[y1] <- sapply(y1, function(y) ((runif(1) > x$p[h]) + 0L))

      rm(y1)
    }

    cohort.freqs <- data.frame(Sum = tapply(haplo.frame$Mum.Allele, haplo.frame$cohort, sum) + tapply(haplo.frame$Dad.Allele, haplo.frame$cohort, sum),
                               Count = tapply(haplo.frame$cohort, haplo.frame$cohort, length))

    cohort.freqs$p <- 1 - (cohort.freqs$Sum/(cohort.freqs$Count * 2))

    cohort.freqs$Sum <- NULL
    cohort.freqs$Count <- NULL

    #~~ sample the rest

    for(h in (n_founder_cohorts+1):nrow(x)){

      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))
      haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")], 1, function(y) y[((runif(1) > 0.5) + 1L)])

      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))
      haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")], 1, function(y) y[((runif(1) > 0.5) + 1L)])


      cohort.freqs$p[h] <- 1 - (sum(haplo.frame$Mum.Allele[y1]) + sum(haplo.frame$Dad.Allele[y2]))/(length(y1) + length(y2))

      if(is.na(cohort.freqs$p[h])){
        stop(paste("Cohort frequency can't be estimated. Problem simulation", simulation, "generation", h))
      }

      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$MOTHER))
      haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) ((runif(1) > cohort.freqs$p[h]) + 0L))

      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$FATHER))
      haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) ((runif(1) > cohort.freqs$p[h]) + 0L))

      rm(y1, y2, y3, y4)
    }

    head(haplo.frame)

    haplo.frame$MOTHER <- NULL
    haplo.frame$FATHER <- NULL

    haplo.frame$Simulation <- simulation

    sim.list[[simulation]] <- haplo.frame

    rm(cohort.freqs, haplo.frame, h)

  }

  sim.results <- do.call(rbind, sim.list)
  sim.results$Simulated.Geno <- sim.results$Mum.Allele + sim.results$Dad.Allele
  names(sim.results)[names(sim.results) == "genotype"] <- "True.Geno"

  sim.results

}



