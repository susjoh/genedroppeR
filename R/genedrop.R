#' Conduct a single Gene-Drop Simulation.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id. Missing mothers
#'   should be coded as 0 or NA
#' @param father vector. Paternal IDs corresponding to id. Missing fathers
#'   should be coded as 0 or NA
#' @param cohort vector. Cohort IDs corresponding to id.
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
#' @param verbose logical. Default = TRUE. Output progress of the function.
#' @param return_genotypes logical. Default = TRUE. Output simulated genotypes for
#'   each individual. If F, only returns the allele frequencies.
#' @param keep_only_genotyped_ids logical. Default = TRUE. Keeps only genotyped IDs
#'   to allow a direct comparison of changes in allele frequencies.
#' @import plyr
#' @import reshape2
#' @import kinship2
#' @import magrittr
#' @import svMisc
#' @export


genedrop <- function(id,
                     mother,
                     father,
                     cohort,
                     genotype,
                     nsim,
                     n_founder_cohorts,
                     fix_founders = T,
                     verbose = T,
                     return_genotypes = T,
                     keep_only_genotyped_ids = T,
                     delim = ""){


  require(reshape2)
  require(plyr)
  require(kinship2)
  require(magrittr)
  require(svMisc)

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ FORMAT THE PEDIGREE        #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  ped <- data.frame(ID     = id,
                    MOTHER = mother,
                    FATHER = father)

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA


  #~~ Get rid of parents that aren't in the IDs

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA

  #~~ Check that ID has not been duplicated

  if(any(as.numeric(names(table(table(ped$ID)))) > 1)) stop("Duplicated values in ID column")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ FORMAT THE COHORTS & GENOTYPES #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ deal with cohorts

  if(is.null(cohort)) {

    ped$cohort <- kindepth(ped[,1], ped$FATHER, ped$MOTHER)

  } else {

    ped$cohort <- cohort
  }

  #~~ deal with genotypes

  ped$genotype <- genotype

  #~~ add a founder column to determine how many preceding generations there are.

  ped$founder <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)


  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ CREATE AN OBJECT TO HOLD GENOTYPES   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ melt pedfile to get a unique row for each gamete transfer

  transped <- melt(ped[,c("ID", "MOTHER", "FATHER")], id.vars = "ID")
  transped$variable <- as.character(transped$variable)

  suppressMessages(transped <- join(transped, subset(ped, select = c(ID, cohort, founder))))

  #~~ Redefine columns

  names(transped) <- c("Offspring.ID", "Parent.ID.SEX", "Parent.ID", "Cohort", "Founder")


  #~~ Get the founder frequency information

  x <- locus_summary(ped$ID, ped$MOTHER, ped$FATHER, ped$cohort, ped$genotype, delim = delim)
  x <- na.omit(x)

  ped <- subset(ped, cohort %in% x$cohort)

  allele.ids <- gsub("Freq_", "", grep("Freq_", names(x), value = T))

  if(length(allele.ids) == 1) stop("Locus is monomorphic. Genedrop simulation terminated")

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ If founders are fixed then create data frame with their information   #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if(fix_founders == T){

    temp <- subset(ped, !is.na(ped$genotype))
    temp$AlleleMOTHER <- sapply(temp$genotype, function (bar) strsplit(bar, split = delim)[[1]][1])
    temp$AlleleFATHER <- sapply(temp$genotype, function (bar) strsplit(bar, split = delim)[[1]][2])

  } else {

    temp <- data.frame(ID = NA)

  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  #~~ Sample the genotypes in the founder cohorts                           #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  if(keep_only_genotyped_ids) return_genotypes <- T

  sim.results <- list()

  for(simulation in 1:nsim){

    if(verbose) progress(simulation, max.value = nsim)


    #~~ Create a list to sample haplotypes

    haplo.list <- list()

    haplo.list[1:length(unique(transped$Offspring.ID))] <- list(list(MOTHER = NA, FATHER = NA))

    names(haplo.list) <- unique(as.character(transped$Offspring.ID))


    for(h in 1:n_founder_cohorts){

      if(verbose == TRUE & nsim == 1) message(paste("Sampling founder genotypes for cohort", x$cohort[h]))

      for(i in which(transped$Cohort == x$cohort[h])){

        allele.founder.freqs <- unlist(x[h, grep("Freq_", names(x))])


        if(transped$Parent.ID.SEX[i] == "MOTHER" & !transped$Offspring.ID[i] %in% temp$ID) {

          haplo.list[as.character(transped$Offspring.ID[i])][[1]]$MOTHER <- sample(allele.ids,
                                                                                   size = 1,
                                                                                   replace = T,
                                                                                   prob = allele.founder.freqs)

        }

        if(transped$Parent.ID.SEX[i] == "FATHER" & !transped$Offspring.ID[i] %in% temp$ID) {

          haplo.list[as.character(transped$Offspring.ID[i])][[1]]$FATHER <- sample(allele.ids,
                                                                                   size = 1,
                                                                                   replace = T,
                                                                                   prob = allele.founder.freqs)

        }


        if(transped$Parent.ID.SEX[i] == "MOTHER" & transped$Offspring.ID[i] %in% temp$ID) {

          temp$AlleleMOTHER[which(temp$ID == transped$Offspring.ID[i])]

          haplo.list[as.character(transped$Offspring.ID[i])][[1]]$MOTHER <- temp$AlleleMOTHER[which(temp$ID == transped$Offspring.ID[i])]


        }

        if(transped$Parent.ID.SEX[i] == "FATHER" & transped$Offspring.ID[i] %in% temp$ID) {

          haplo.list[as.character(transped$Offspring.ID[i])][[1]]$FATHER <- temp$AlleleFATHER[which(temp$ID == transped$Offspring.ID[i])]

        }

      }

    }


    #~~ Create holders for the simulated frequencies

    cohort.freqs <- NULL

    for(h in 1:n_founder_cohorts){

      y <- haplo.list[as.character(ped$ID[which(ped$cohort %in% x$cohort[h])])] %>% unlist
      y <- data.frame(table(y)/length(y))
      y <- cbind(cohort = x$cohort[h], y)
      suppressMessages(y <- dcast(y, cohort ~ y))

      if(ncol(y) < length(allele.ids) + 1) eval(parse(text = paste0("y$", allele.ids[which(!allele.ids %in% names(y))], " <- 0")))

      cohort.freqs <- rbind(cohort.freqs, y)

      rm(y)

    }

    #~~ Sample the non-founder haplotypes by cohort. This loops through cohorts sequentially as
    #   parental haplotypes must exist before sampling. Sample the founders from the same cohort

    for(i in x$cohort[(n_founder_cohorts+1):nrow(x)]){

      if(verbose == TRUE & nsim == 1) message(paste("Simulating genotypes for cohort", i))

      #~~ Run the IDs with known parents.

      for(j in which(transped$Cohort == i & !is.na(transped$Parent.ID))){

        #~~ Sample one of the parental alleles at random

        haplo.list[as.character(transped$Offspring.ID[j])][[1]][transped$Parent.ID.SEX[j]] <- haplo.list[as.character(transped$Parent.ID[j])][[1]][(runif(1) < 0.5) + 1L]

      }


      y <- haplo.list[as.character(ped$ID[which(ped$cohort %in% i)])] %>% unlist
      y <- na.omit(y)
      y <- factor(y, levels = allele.ids)
      allele.founder.freqs <- table(y)/length(y)
      rm(y)

      #~~ Now sample the founders from the rest of the cohort

      for(k in which(transped$Cohort == i & is.na(transped$Parent.ID))){

        if(transped$Parent.ID.SEX[k] == "MOTHER") {

          haplo.list[as.character(transped$Offspring.ID[k])][[1]]$MOTHER <- sample(allele.ids,
                                                                                   size = 1,
                                                                                   replace = T,
                                                                                   prob = allele.founder.freqs)

        }

        if(transped$Parent.ID.SEX[k] == "FATHER") {

          haplo.list[as.character(transped$Offspring.ID[k])][[1]]$FATHER <- sample(allele.ids,
                                                                                   size = 1,
                                                                                   replace = T,
                                                                                   prob = allele.founder.freqs)

        }

      }



      y <- haplo.list[as.character(ped$ID[which(ped$cohort %in% i)])] %>% unlist
      y <- data.frame(table(y)/length(y))
      y <- cbind(cohort = i, y)
      suppressMessages(y <- dcast(y, cohort ~ y))

      if(ncol(y) < length(allele.ids) + 1) eval(parse(text = paste0("y$", allele.ids[which(!allele.ids %in% names(y))], " <- 0")))


      cohort.freqs <- rbind(cohort.freqs, y)

      rm(y)

    }

    cohort.freqs$buildpop <- "simulated"
    cohort.freqs$buildpop[1:n_founder_cohorts] <- "sampled"

    cohort.freqs <- melt(cohort.freqs, id.vars = c("cohort", "buildpop"))
    names(cohort.freqs)[3:4] <- c("Allele", "Frequency")


    if(return_genotypes == T){

      if(verbose & nsim == 1) message("Formatting genotype data...")

      genotype.list <- lapply(haplo.list, function(foo) data.frame(Allele1 = foo[1][[1]],
                                                                   Allele2 = foo[2][[1]]))

      geno.table <- do.call(rbind, genotype.list)

      geno.table$ID <- row.names(geno.table)

      suppressMessages(geno.table <- join(geno.table, subset(ped, select = c(ID, cohort))))

      newlist <- list(cohort.freqs, geno.table)

      names(newlist) <- c("cohort.freqs", "geno.table")


    } else {

      newlist <- list(cohort.freqs)

      names(newlist) <- c("cohort.freqs")


    }

    if(verbose & nsim == 1) message("...done.")

    sim.results[[simulation]] <- newlist
    sim.results[[simulation]]$cohort.freqs$Iteration <- simulation
    if(return_genotypes) sim.results[[simulation]]$geno.table$Iteration <- simulation

  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # FORMAT THE RESULTS                          #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

  #~~ extract cohort genotypes

  if(keep_only_genotyped_ids){

    #~~ get the genotypes out

    sim.genotypes <- NULL

    sim.genotypes <- lapply(sim.results, function(x) x$geno.table)

    sim.genotypes <- do.call(rbind, sim.genotypes)

    #~~ Which IDs are genotyped?

    sim.genotypes <- subset(sim.genotypes, ID %in% subset(ped, !is.na(genotype))$ID) %>% droplevels
    sim.genotypes$Allele1 <- as.character(sim.genotypes$Allele1)
    sim.genotypes$Allele2 <- as.character(sim.genotypes$Allele2)


    #~~ Get the allele frequencies per year.

    sim.genotypes <- melt(sim.genotypes, id.vars = c("ID", "cohort", "Iteration"))

    sim.alleles.summary   <- data.frame(table(sim.genotypes$Iteration, sim.genotypes$cohort, sim.genotypes$value))
    sim.genotypes.summary <- data.frame(table(sim.genotypes$Iteration, sim.genotypes$cohort))

    names(sim.alleles.summary) <- c("Iteration", "cohort", "Allele", "Freq")
    names(sim.genotypes.summary) <- c("Iteration", "cohort", "Total")

    suppressMessages(sim.alleles.summary <- join(sim.alleles.summary, sim.genotypes.summary))
    sim.alleles.summary$Frequency <- sim.alleles.summary$Freq/sim.alleles.summary$Total

    sim.alleles.summary$Freq <- NULL
    sim.alleles.summary$Total <- NULL

    sim.alleles.summary$buildpop <- "simulated"
    sim.alleles.summary$buildpop[which(sim.alleles.summary$cohort %in% x$cohort[1:n_founder_cohorts])] <- "sampled"

    sim.alleles.summary <- sim.alleles.summary[,c("cohort", "buildpop", "Allele", "Frequency", "Iteration")]

    rm(sim.genotypes, sim.genotypes.summary)

  } else {

    sim.alleles.summary <- lapply(sim.results, function(x) x$cohort.freqs)
    sim.alleles.summary <- do.call(rbind, sim.alleles.summary)

  }

  sim.alleles.summary$cohort <- as.numeric(as.character(sim.alleles.summary$cohort))
  sim.alleles.summary$Iteration <- as.numeric(as.character(sim.alleles.summary$Iteration))


  #~~ Extract the true results

  genotypeSummaryFreqs <- function(genotypeSummaryObject){

    require(reshape2)
    x <- genotypeSummaryObject[,c("cohort", grep("Freq", names(genotypeSummaryObject), value = T))]
    x <- melt(x, id.vars = "cohort")
    names(x)[2:3] <- c("Allele", "Frequency")

    x$Allele <- gsub("Freq_", "", x$Allele)

    x
  }

  true.alleles.summary <- x %>% genotypeSummaryFreqs
  true.alleles.summary$Iteration = 0

  true.alleles.summary$buildpop <- "simulated"
  true.alleles.summary$buildpop[which(true.alleles.summary$cohort %in% x$cohort[1:n_founder_cohorts])] <- "sampled"

  true.alleles.summary <- true.alleles.summary[,c("cohort", "buildpop", "Allele", "Frequency", "Iteration")]

  #~~ create a genotype table if required

  if(return_genotypes){

    sim.genotypes <- lapply(sim.results, function(x) x$geno.table)
    sim.genotypes <- do.call(rbind, sim.genotypes)

  }

  #~~ return results

  if(return_genotypes){

    x <- list(sim.alleles.summary, true.alleles.summary, sim.genotypes)
    names(x) <- c("cohort.freqs", "observed.freqs", "geno.table")

  } else {

    x <- list(sim.alleles.summary, true.alleles.summary)
    names(x) <- c("cohort.freqs", "observed.freqs")


  }

  return(x)



}
