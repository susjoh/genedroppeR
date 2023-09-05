#' Plot a summary graph of Gene-Drop Simulations.
#'
#' Plot a summary graph of Gene-Drop Simulations.
#'
#' @param genedrop_object_summary Gene-Drop summary object from the function `summary_genedrop()`
#' @param sim_alpha alpha (transparency) value for plotted points
#' @param obs_line_col line colour to use for the observed data.
#' @import ggplot2
#' @export



plot_genedrop_results <- function(genedrop_object_summary,
                                  sim_alpha = 0.2,
                                  obs_line_col = "red"){

  genedrop_object_summary$simulated_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Simulation))
  genedrop_object_summary$simulated_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Cohort))

  genedrop_object_summary$observed_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Simulation))
  genedrop_object_summary$observed_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Cohort))


  if(!"Allele" %in% names(genedrop_object_summary$simulated_frequencies)){
    genedrop_object_summary$simulated_frequencies$Allele <- "p"
    genedrop_object_summary$observed_frequencies$Allele  <- "p"
  }


  ggplot(genedrop_object_summary$simulated_frequencies, aes(Cohort, p, group = Simulation)) +
    geom_line(alpha = sim_alpha) +
    geom_line(data = genedrop_object_summary$observed_frequencies, aes(Cohort, p), col = obs_line_col) +
    ggtitle(paste0("Allele Frequency Changes: Nsim = ", max(genedrop_object_summary$simulated_frequencies$Simulation))) +
    facet_wrap(~Allele)

}

#' Plot a histogram of Gene-Drop Simulation linear regression slopes
#'
#' Plot a histogram of Gene-Drop Simulation linear regression slopes and return
#' true distribution values. #'
#' @param genedrop_object_summary Gene-Drop summary object from the function
#'   `summary_genedrop()`
#' @param n_founder_cohorts integer. The number of cohorts at the top of the
#'   pedigree that will sample from the true allele frequences (these are
#'   defined as "sampled"). All cohorts following these ones are "simulated" and
#'   are used for comparisons of changes in allele frequency.
#' @param remove_founders logical. Default = TRUE. Remove the founders from the
#'   data before calculating the linear regression slopes.
#' @param method Default = "lm". Currently the only option.
#' @param obs_line_col line colour to use for the observed data.
#' @import ggplot2
#' @export


plot_genedrop_lm_slopes <- function(genedrop_object_summary,
                                    n_founder_cohorts = NULL,
                                    remove_founders = T,
                                    method = "lm",
                                    obs_line_col = "red"){

  sim.slopes <- NULL

  genedrop_object_summary$simulated_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Simulation))
  genedrop_object_summary$simulated_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Cohort))

  genedrop_object_summary$observed_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Simulation))
  genedrop_object_summary$observed_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Cohort))

  if(remove_founders){
    if(length(n_founder_cohorts) == 1){

      genedrop_object_summary$simulated_frequencies <-
        subset(genedrop_object_summary$simulated_frequencies,
               !Cohort %in% unique(sort(genedrop_object_summary$simulated_frequencies$Cohort))[1:n_founder_cohorts])

      genedrop_object_summary$observed_frequencies <-
        subset(genedrop_object_summary$observed_frequencies,
               !Cohort %in% unique(sort(genedrop_object_summary$observed_frequencies$Cohort))[1:n_founder_cohorts])

    }
  }

  if("Allele" %in% names(genedrop_object_summary$simulated_frequencies)){

    Allele = unique(genedrop_object_summary$simulated_frequencies$Allele)

  } else {

    Allele = "p"

  }

  for(i in 1:max(genedrop_object_summary$simulated_frequencies$Simulation)){

    if(i %in% seq(1, max(genedrop_object_summary$simulated_frequencies$Simulation), 100)) print(paste("Calculating slope", i, "in", max(genedrop_object_summary$simulated_frequencies$Simulation)))

    for(j in Allele){

      x <- subset(genedrop_object_summary$simulated_frequencies, Simulation == i & Allele == j)
      x1 <- lm(p ~ Cohort, data = x)$coefficients[[2]]

      sim.slopes <- rbind(sim.slopes,
                          data.frame(Iteration = i,
                                     Allele = j,
                                     Slope = x1))

      rm(x1)

    }

  }

  true.slopes <- NULL

  for(i in 0){

    for(j in Allele){

      x <- subset(genedrop_object_summary$observed_frequencies, Simulation == i & Allele == j)
      x1 <- lm(p ~ Cohort, data = x)$coefficients[[2]]

      true.slopes <- rbind(true.slopes,
                           data.frame(Iteration = i,
                                      Allele = j,
                                      Slope = x1))

      rm(x1)

    }

  }


  print(ggplot(sim.slopes, aes(Slope)) +
          geom_histogram() +
          facet_wrap(~Allele) +
          geom_vline(data = true.slopes, aes(xintercept = Slope), col = obs_line_col) +
          ggtitle(paste0("Distribution of Regression Slopes: Nsim = ", max(genedrop_object_summary$simulated_frequencies$Simulation))))

  true.slopes$Slopes.Lower <- NA
  true.slopes$Slopes.Higher <- NA

  for(i in 1:nrow(true.slopes)){

    true.slopes$Slopes.Lower [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                   sim.slopes$Slope < true.slopes$Slope[i]))

    true.slopes$Slopes.Higher [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                    sim.slopes$Slope > true.slopes$Slope[i]))

  }

  true.slopes

}


#' Plot a histogram of Gene-Drop Simulation cumulative change values.
#'
#' Plot a histogram of Gene-Drop Simulation cumulated change values and return
#' true distribution values. This is the cumulative change in allele frequencies
#' over time and may be a signature of balancing selection.
#'
#' @param genedrop_object_summary Gene-Drop summary object from the function
#'   `summary_genedrop()`
#' @param n_founder_cohorts integer. The number of cohorts at the top of the
#'   pedigree that will sample from the true allele frequences (these are
#'   defined as "sampled"). All cohorts following these ones are "simulated" and
#'   are used for comparisons of changes in allele frequency.
#' @param remove_founders logical. Default = TRUE. Remove the founders from the
#'   data before calculating the linear regression slopes.
#' @param method Default = "lm". Currently the only option.
#' @param obs_line_col line colour to use for the observed data.
#' @import ggplot2
#' @export


plot_genedrop_cumulative_change <- function(genedrop_object_summary,
                                            n_founder_cohorts = NULL,
                                            remove_founders = T,
                                            obs_line_col = "red"){


  genedrop_object_summary$simulated_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Simulation))
  genedrop_object_summary$simulated_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$simulated_frequencies$Cohort))

  genedrop_object_summary$observed_frequencies$Simulation <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Simulation))
  genedrop_object_summary$observed_frequencies$Cohort <- as.numeric(as.character(genedrop_object_summary$observed_frequencies$Cohort))

  if(remove_founders){
    if(length(n_founder_cohorts) == 1){

      genedrop_object_summary$simulated_frequencies <-
        subset(genedrop_object_summary$simulated_frequencies,
               !Cohort %in% unique(sort(genedrop_object_summary$simulated_frequencies$Cohort))[1:n_founder_cohorts])

      genedrop_object_summary$observed_frequencies <-
        subset(genedrop_object_summary$observed_frequencies,
               !Cohort %in% unique(sort(genedrop_object_summary$observed_frequencies$Cohort))[1:n_founder_cohorts])

    }
  }

  if("Allele" %in% names(genedrop_object_summary$simulated_frequencies)){

    Allele = unique(genedrop_object_summary$simulated_frequencies$Allele)

  } else {

    Allele = "p"
    genedrop_object_summary$simulated_frequencies$Allele <- "p"
    genedrop_object_summary$observed_frequencies$Allele <- "p"

  }

  Allele <- as.character(Allele)

  cumu.func <- function(x){
    x <- diff(x)
    x <- ifelse(x < 0, x * -1, x)
    sum(x, na.rm = F)
  }

  sim.changes <- tapply(genedrop_object_summary$simulated_frequencies$p,
                        list(genedrop_object_summary$simulated_frequencies$Simulation,
                             genedrop_object_summary$simulated_frequencies$Allele),
                        cumu.func)

  sim.changes <- melt(sim.changes)

  true.changes <- tapply(genedrop_object_summary$observed_frequencies$p,
                        list(genedrop_object_summary$observed_frequencies$Simulation,
                             genedrop_object_summary$observed_frequencies$Allele),
                        cumu.func)

  true.changes <- melt(true.changes)



  print(ggplot(sim.changes, aes(value)) +
          geom_histogram() +
          facet_wrap(~Var2) +
          geom_vline(data = true.changes, aes(xintercept = value), col = obs_line_col) +
          ggtitle(paste0("Distribution of Cumulative Change: Nsim = ", max(genedrop_object_summary$simulated_frequencies$Simulation))))

  true.changes$Cumulative.Change.Lower <- NA
  true.changes$Cumulative.Change.Higher <- NA

  for(i in 1:nrow(true.changes)){

    true.changes$Cumulative.Change.Lower [i] <- length(which(sim.changes$Var2 == true.changes$Var2[i] &
                                                               sim.changes$value < true.changes$value[i]))

    true.changes$Cumulative.Change.Higher [i] <- length(which(sim.changes$Var2 == true.changes$Var2[i] &
                                                                sim.changes$value > true.changes$value[i]))

  }

  true.changes

}


