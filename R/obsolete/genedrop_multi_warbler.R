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
#' @param verbose logical. Output the progress of the run.
#' @param interval int. Default 100. Output progress every 100 simulations.




genedrop_multi_warbler <- function(id,
                                   mother,
                                   father,
                                   cohort = NULL,
                                   genotype,
                                   genotype_delim = '',
                                   nsim,
                                   n_founder_cohorts = 1,
                                   fix_founders = T,
                                   verbose = T,
                                   interval = 100){

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Format the data           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ Check that there are no duplicate IDs.

  if (any(as.numeric(names(table(table(id)))) > 1)){
    stop ("Duplicated values in id")
  }

  #~~ If cohort is provided, then check there are no NA's.

  if (!is.null(cohort) & any(is.na(cohort))){
    message("NAs present in cohort - these individuals will be removed.")
  }

  #~~ If genotype is provided, then check the locus is not monomorphic and is biallelic.

  if (!is.null(genotype) & length(table(genotype)) == 1){
    stop ("Locus is monomorphic")
  }

  if(max(sapply(y$genotype, function(foo) length(strsplit(foo, split = genotype_delim, fixed = T)[[1]]))) > 2){
    stop ("There are more than two alleles in some IDs - check genotype_delim.")
  }


  #~~ Make a ped object.

  ped <- data.frame(ID     = id,
                    MOTHER = mother,
                    FATHER = father)

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA

  #~~ Add the genotype information.

  if (!is.null(genotype)) ped$genotype <- as.character(genotype)

  #~~ Add cohort information to ped.

  if (!is.null(cohort)) {

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

    if (length(badids) > 0){
      print(x[badids,])
      warning ("Some parents and offspring are in the same cohort.
            See output for problem lines.
            Model will run despite this.")
    }

    rm(x, x1, badids)

    # deal with fact that two IDs could be in the same cohort

    temp <- NULL

    for(h in 1:max(ped$cohort)){

      temp2 <- subset(ped, cohort == h)
      temp2$cohort2 <- kindepth(temp2$ID, temp2$MOTHER, temp2$FATHER)
      temp <- rbind(temp, temp2)
      rm(temp2)
    }

    ped <- temp

    rm(temp)
  } else {

    ped$cohort <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)
    ped$cohort2 <- 0


  }

  #~~ Get rid of parents that aren't in the IDs.

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA


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
  head(y)

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

  sim.list <- list()

  for(simulation in 1:nsim){

    if (verbose){
      if (simulation %in% seq(1, nsim, interval)){
        message(paste0("Running simulation ", simulation, " of ", nsim, "."))
      }
    }

    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    #~~ Sample the founders

    for(h in 1:n_founder_cohorts){

      for(j in sort(unique(subset(haplo.frame, cohort == h)$cohort2))){

        y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER) & !haplo.frame$ID %in% fixed_founders & haplo.frame$cohort2 == j)

        haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                            1,
                                            function(y) y[((runif (1) > 0.5) + 1L)])

        y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER) & !haplo.frame$ID %in% fixed_founders & haplo.frame$cohort2 == j)

        haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                            1,
                                            function(y) y[((runif (1) > 0.5) + 1L)])

        y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Mum.Allele) & haplo.frame$cohort2 == j)

        if(length(y3) > 0){
          haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) sample(x.allele, size = 1, prob = x[h, x.allele]))
        }

        y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Dad.Allele) & haplo.frame$cohort2 == j)

        if(length(y4) > 0){
          haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) sample(x.allele, size = 1, prob = x[h, x.allele]))
        }
        rm(y1, y2, y3, y4)
      }
    }

    #~~ sample the rest

    for(h in (n_founder_cohorts+1):nrow(x)){

      for(j in sort(unique(subset(haplo.frame, cohort == h)$cohort2))){

        y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER) & haplo.frame$cohort2 == j)

        haplo.frame$Mum.Allele[y1] <- apply(haplo.frame[haplo.frame$MOTHER[y1],c("Mum.Allele", "Dad.Allele")],
                                            1,
                                            function(y) y[((runif (1) > 0.5) + 1L)])

        y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER)  & haplo.frame$cohort2 == j)

        haplo.frame$Dad.Allele[y2] <- apply(haplo.frame[haplo.frame$FATHER[y2],c("Mum.Allele", "Dad.Allele")],
                                            1,
                                            function(y) y[((runif (1) > 0.5) + 1L)])

        rm(y1, y2)
      }
      #~~ Get allele frequencies

      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$MOTHER))

      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$FATHER))


      temp.freq <- data.frame(table(c(haplo.frame$Mum.Allele[y1], haplo.frame$Dad.Allele[y2])))
      temp.freq$Freq <- temp.freq$Freq / sum(temp.freq$Freq)
      temp.freq$Var1 <- as.character(temp.freq$Var1)

      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$MOTHER))

      haplo.frame$Mum.Allele[y3] <- sapply(y3, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$FATHER))

      haplo.frame$Dad.Allele[y4] <- sapply(y4, function(y) sample(temp.freq$Var1, size = 1, prob = temp.freq$Freq))

      rm(y3, y4)
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

  names(sim.results)[names(sim.results) == "genotype"] <- "True.Geno"

  sim.results

}


