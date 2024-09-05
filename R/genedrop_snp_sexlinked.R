#' `genedrop_snp_sexlinked()`: Conduct a genedrop simulation for a sex-linked
#' biallelic locus.
#'
#' This function conducts a genedrop simulation for a single, sex-linked
#' bi-allelic locus (e.g. a, X- or Z-linked SNP). For autosomal loci, use
#' `genedrop_snp()`. Before running this function, users should first summarise
#' and visualise their data using `summary_cohort()` to determine an appropriate
#' value for `n_founder_cohorts`. This function will return an object that
#' contains the cohort allele frequences in the observed and simulated datasets.
#' Overall results of directional and balancing selection can be observed using
#' `summary()`. For more detail on specifying model parameters, please consult
#' the tutorial at https://github.com/susjoh/genedroppeR.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param sex vector. Sexes corresponding to id. 1 is the homogametic sex (e.g.
#'   XY, ZW) and 2 is the heterogametic sex (e.g. XX, ZZ)
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
#' @param remove_founders Default = TRUE. If TRUE, then the founder cohorts will
#'   be removed from calculations of directional and cumulative change.
#' @param return_full_results Default = NULL. This will also output tables of
#'   all individually simulated genotypes.
#' @param verbose logical. Default = TRUE. Output the progress of the run.
#' @param interval integer. Default 100. Output progress every 100 simulations.
#' @examples
#' data(unicorn)
#' sub_unicorn <- subset(unicorn, sex %in% c(1,2))
#' genedrop_obj <- genedrop_snp_sex(id = sub_unicorn$id,
#'                                  mother = sub_unicorn$mother,
#'                                  father = sub_unicorn$father,
#'                                  cohort = sub_unicorn$cohort,
#'                                  genotype = sub_unicorn$Xlinked,
#'                                  sex = sub_unicorn$sex,
#'                                  nsim = 100,
#'                                  n_founder_cohorts = 4,
#'                                  fix_founders = TRUE,
#'                                  verbose = TRUE,
#'                                  interval = 10)
#' summary_genedrop(genedrop_obj)
#' plot_genedrop(genedrop_obj)
#' @export

genedrop_snp_sex <- function(id,
                         mother,
                         father,
                         cohort = NULL,
                         genotype,
                         sex,
                         nsim,
                         n_founder_cohorts = 1,
                         fix_founders = TRUE,
                         verbose = TRUE,
                         interval = 100,
                         resample_offspring = FALSE,
                         remove_founders = TRUE,
                         return_full_results = NULL){

  Hom.Parent.Allele = Het.Parent.Allele = Cohort = Simulation = p = NULL

  # Check the data and obtain ped object

  ped_check <- check_data(id, mother, father, cohort, genotype, sex)
  ped <- ped_check$ped
  sex_system <- ped_check$sex_system

  rm(id, mother, father, cohort, genotype, sex, ped_check)

  #~~ Recode to Het and Hom parent.

  if(sex_system == "XY"){
    ped$HET_Parent = ped$FATHER
    ped$HOM_Parent = ped$MOTHER
  } else {
    ped$HET_Parent = ped$MOTHER
    ped$HOM_Parent = ped$FATHER
  }

  # deal with heterogametic sex errors

  heterrors <- which(ped$genotype == 1 & ped$sex == 1)

  if(length(heterrors) > 0){
    message(paste0("Removed ", length(heterrors), " heterozygous genotypes in ", sex_system, " ids (", round(((length(heterrors)/length(which(ped$sex == 1)))*100), 3), "%)."))
    ped$genotype[heterrors] <- NA
  }

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

    #~~ index the pedigree

    row.names(ped) <- ped$ID

    #~~ Create columns for parentally inherited alleles and add alleles for all
    #   individuals that are in the founder cohorts.

    ped$Hom.Parent.Allele <- NA
    ped$Het.Parent.Allele <- NA

    if(fix_founders){

      # Generate alleles
      ped$Hom.Parent.Allele <- ifelse(ped$genotype %in% 0:1, 0, ifelse(is.na(ped$genotype), NA, 1))
      ped$Het.Parent.Allele <- ifelse(ped$genotype %in% 0  , 0, ifelse(is.na(ped$genotype), NA, 1))

      # Blank out anything that is not in the founder cohorts
      ped$Hom.Parent.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA
      ped$Het.Parent.Allele[which(ped$cohort %in% x$cohort[(n_founder_cohorts+1):nrow(x)])] <- NA

      # Blank out founders that have to inherit an allele from a parent
      ped$Hom.Parent.Allele[which(!is.na(ped$HOM_Parent))] <- NA
      ped$Het.Parent.Allele[which(!is.na(ped$HET_Parent))] <- NA

    }

    # Convert to character

    ped$ID <- as.character(ped$ID)
    ped$HOM_Parent <- as.character(ped$HOM_Parent)
    ped$HET_Parent <- as.character(ped$HET_Parent)

    #~~ Create a data frame with space for the results

    haplo.frame <- ped

    #~~ Sample the founders

    for(h in 1:n_founder_cohorts){

      # HOM_Parent is known
      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$HOM_Parent))

      if(length(y1) > 0)  haplo.frame$Hom.Parent.Allele[y1] <- apply(haplo.frame[haplo.frame$HOM_Parent[y1],c("Hom.Parent.Allele", "Het.Parent.Allele")],
                                                                     1,
                                                                     function(y) y[((runif (1) > 0.5) + 1L)])

      # HET_Parent is known
      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$HET_Parent) & haplo.frame$sex == 2)

      if(length(y2) > 0)  haplo.frame$Het.Parent.Allele[y2] <- haplo.frame[haplo.frame$HET_Parent[y2],"Hom.Parent.Allele"]

      # HOM_Parent allele is not known - sample from cohort frequency
      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Hom.Parent.Allele))

      if(length(y3) > 0){
        haplo.frame$Hom.Parent.Allele[y3] <- sapply(y3, function(y) ((runif (1) > x$p[h]) + 0L))
      }

      # HET_Parent allele is not known - sample from cohort frequency
      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$Het.Parent.Allele) & haplo.frame$sex == 2)

      if(length(y4) > 0){
        haplo.frame$Het.Parent.Allele[y4] <- sapply(y4, function(y) ((runif (1) > x$p[h]) + 0L))
      }

      # Which heterogametics have a value for their heterogametic parent? Make them have the same as their homogametic parent

      haplo.frame$Het.Parent.Allele[which(haplo.frame$sex == 1)] <- NA

      rm(y1, y2, y3, y4)

    }

    #~~ Calculate the cohort frequencies

    cohort.freqs <- haplo.frame %>% group_by(cohort) %>% summarise(Sum = sum(Hom.Parent.Allele, Het.Parent.Allele, na.rm = T),
                                                                   Count = length(na.omit(c(Hom.Parent.Allele, Het.Parent.Allele))))

    cohort.freqs$p <- 1 - (cohort.freqs$Sum / (cohort.freqs$Count * 2))

    cohort.freqs$Sum <- NULL
    cohort.freqs$Count <- NULL

    #~~ sample the rest

    for(h in (n_founder_cohorts+1):nrow(x)){

      # HOM_Parent is known
      y1 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$HOM_Parent))

      if(length(y1) > 0)  haplo.frame$Hom.Parent.Allele[y1] <- apply(haplo.frame[haplo.frame$HOM_Parent[y1],c("Hom.Parent.Allele", "Het.Parent.Allele")],
                                                                     1,
                                                                     function(y) y[((runif (1) > 0.5) + 1L)])

      # HET_Parent is known
      y2 <- which(haplo.frame$cohort == x$cohort[h] & !is.na(haplo.frame$HET_Parent) & haplo.frame$sex == 2)

      if(length(y2) > 0)  haplo.frame$Het.Parent.Allele[y2] <- haplo.frame[haplo.frame$HET_Parent[y2],"Hom.Parent.Allele"]


      # Estimate the allele frequency
      cohort.freqs$p[h] <- 1 - (sum(haplo.frame$Hom.Parent.Allele[y1]) + sum(haplo.frame$Het.Parent.Allele[y2])) / (length(y1) + length(y2))

      if (is.na(cohort.freqs$p[h])){
        stop (paste("Cohort frequency can't be estimated. Problem simulation", simulation, "generation", h))
      }

      # Now sample IDs with missing HOM_Parent
      y3 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$HOM_Parent))

      if(length(y3) > 0){
        haplo.frame$Hom.Parent.Allele[y3] <- sapply(y3, function(y) ((runif (1) > cohort.freqs$p[h]) + 0L))
      }

      # Now sample IDs with missing HET_Parent
      y4 <- which(haplo.frame$cohort == x$cohort[h] & is.na(haplo.frame$HET_Parent) & haplo.frame$sex == 2)

      if(length(y4) > 0){
        haplo.frame$Het.Parent.Allele[y4] <- sapply(y4, function(y) ((runif (1) > cohort.freqs$p[h]) + 0L))
      }
      rm(y1, y2, y3, y4)

    }

    haplo.frame$HOM_Parent <- NULL
    haplo.frame$HET_Parent <- NULL

    haplo.frame$Simulation <- simulation

    sim.list[[simulation]] <- haplo.frame

    rm(cohort.freqs, haplo.frame, h)

  }

  sim.results <- bind_rows(sim.list)

  sim.results$Simulated.Geno <- sim.results$Hom.Parent.Allele + sim.results$Hom.Parent.Allele

  names(sim.results)[names(sim.results) == "genotype"] <- "True.Geno"

  if(!is.null(return_full_results)) return_full_results <- sim.results

  genedrop_obj <- process_genedrop(sim.results)

  # Calculate selection

  genedrop_obj$simulated_frequencies$Simulation <- as.numeric(as.character(genedrop_obj$simulated_frequencies$Simulation))
  genedrop_obj$simulated_frequencies$Cohort <- as.numeric(as.character(genedrop_obj$simulated_frequencies$Cohort))

  genedrop_obj$observed_frequencies$Simulation <- as.numeric(as.character(genedrop_obj$observed_frequencies$Simulation))
  genedrop_obj$observed_frequencies$Cohort <- as.numeric(as.character(genedrop_obj$observed_frequencies$Cohort))

  sim_freq_hold <- genedrop_obj$simulated_frequencies
  obs_freq_hold <- genedrop_obj$observed_frequencies

  if(remove_founders){
    if(length(n_founder_cohorts) == 1){

      genedrop_obj$simulated_frequencies <-
        filter(genedrop_obj$simulated_frequencies,
               !Cohort %in% unique(sort(genedrop_obj$simulated_frequencies$Cohort))[1:n_founder_cohorts])

      genedrop_obj$observed_frequencies <-
        filter(genedrop_obj$observed_frequencies,
               !Cohort %in% unique(sort(genedrop_obj$observed_frequencies$Cohort))[1:n_founder_cohorts])

    }
  }

  if("Allele" %in% names(genedrop_obj$simulated_frequencies)){

    Allele = unique(genedrop_obj$simulated_frequencies$Allele)

  } else {

    sim_freq_hold$Allele <- "p"
    obs_freq_hold$Allele <- "p"

    genedrop_obj$simulated_frequencies$Allele <- "p"
    genedrop_obj$observed_frequencies$Allele <- "p"
  }

  suppressMessages({

    sim.slopes <- genedrop_obj$simulated_frequencies %>%
      group_by(Simulation, Allele) %>%
      summarise(Estimate = lm(p ~ Cohort)$coefficients[[2]])

    true.slopes <- genedrop_obj$observed_frequencies %>%
      group_by(Simulation, Allele) %>%
      summarise(Estimate = lm(p ~ Cohort)$coefficients[[2]])

    cumu.func <- function(x){
      x <- diff(x)
      x <- ifelse(x < 0, x * -1, x)
      sum(x, na.rm = F)
    }

    sim.changes <- genedrop_obj$simulated_frequencies %>%
      group_by(Simulation, Allele) %>%
      summarise(Estimate = cumu.func(p))

    true.changes <- genedrop_obj$observed_frequencies %>%
      group_by(Simulation, Allele) %>%
      summarise(Estimate = cumu.func(p))

  })

  true.slopes$Estimates.Lower <- NA
  true.slopes$Estimates.Higher <- NA

  for(i in 1:nrow(true.slopes)){

    true.slopes$Estimates.Lower [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                      sim.slopes$Estimate < true.slopes$Estimate[i]))

    true.slopes$Estimates.Higher [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                       sim.slopes$Estimate > true.slopes$Estimate[i]))

  }


  true.changes$Estimates.Lower <- NA
  true.changes$Estimates.Higher <- NA

  for(i in 1:nrow(true.changes)){

    true.changes$Estimates.Lower [i] <- length(which(sim.changes$Allele == true.changes$Allele[i] &
                                                       sim.changes$Estimate < true.changes$Estimate[i]))

    true.changes$Estimates.Higher [i] <- length(which(sim.changes$Allele == true.changes$Allele[i] &
                                                        sim.changes$Estimate > true.changes$Estimate[i]))

  }

  # Table the results

  true.changes$Analysis <- "Cumulative Change"
  true.slopes$Analysis <- "Directional Change"

  restab <- rbind(true.changes, true.slopes)
  restab$Simulation <- nsim
  names(restab)[1] <- "Simulations"
  restab <- restab[,c("Analysis", "Allele", "Estimate", "Estimates.Lower", "Estimates.Higher", "Simulations")]


  # Return results

  sim.results <- list(
    results = restab,
    observed_frequencies  = obs_freq_hold,
    simulated_frequencies = sim_freq_hold,
    full_results          = return_full_results,
    n_founder_cohorts     = n_founder_cohorts,
    remove_founders       = remove_founders,
    fix_founders          = fix_founders,
    resample_offspring    = resample_offspring,
    slopes                = list(true.slopes = true.slopes,
                                 sim.slopes  = sim.slopes),
    cumulative_change = list(true.changes = true.changes,
                             sim.changes  = sim.changes)
  )

  class(sim.results) <- "genedroppeR"

  sim.results

}



