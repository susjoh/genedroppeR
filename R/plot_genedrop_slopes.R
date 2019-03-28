#' Plot Genedrop slopes.
#'
#' Plot slopes from a genedrop analysis
#'
#' @param genedrop_object Output from genedrop() function
#' @import ggplot2
#' @export
#'

plot_genedrop_slopes <- function(genedrop_object){

  require(ggplot2)

  #~~ Make plots

  simulated.frequencies <- genedrop_object$cohort.freqs
  observed.frequencies <- genedrop_object$observed.freqs

  p1 <- ggplot(simulated.frequencies, aes(cohort, Frequency, group = Iteration)) +
    geom_line() +
    facet_wrap(~Allele) +
    geom_line(data = observed.frequencies, aes(cohort, Frequency), col = "red")


  simulated.frequencies <- subset(simulated.frequencies, buildpop == "simulated")
  observed.frequencies  <- subset(observed.frequencies , buildpop == "simulated")

  sim.slopes <- NULL

  for(i in 1:max(simulated.frequencies$Iteration)){

    for(j in unique(simulated.frequencies$Allele)){

      x <- subset(simulated.frequencies, Iteration == i & Allele == j)
      x1 <- lm(Frequency ~ cohort, data = x)$coefficients[[2]]

      sim.slopes <- rbind(sim.slopes,
                          data.frame(Iteration = i,
                                     Allele = j,
                                     Slope = x1))

      rm(x1)

    }

  }

  true.slopes <- NULL

  for(i in 0){

    for(j in unique(observed.frequencies$Allele)){

      x <- subset(observed.frequencies, Iteration == i & Allele == j)
      x1 <- lm(Frequency ~ cohort, data = x)$coefficients[[2]]

      true.slopes <- rbind(true.slopes,
                           data.frame(Iteration = i,
                                      Allele = j,
                                      Slope = x1))

      rm(x1)

    }

  }


  p2 <- ggplot(sim.slopes, aes(Slope)) +
    geom_histogram(binwidth = 0.001) +
    facet_wrap(~Allele) +
    geom_vline(data = true.slopes, aes(xintercept = Slope), col = "red")

  true.slopes$Slopes.Lower <- NA
  true.slopes$Slopes.Higher <- NA

  for(i in 1:nrow(true.slopes)){

    true.slopes$Slopes.Lower [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                   sim.slopes$Slope < true.slopes$Slope[i]))

    true.slopes$Slopes.Higher [i] <- length(which(sim.slopes$Allele == true.slopes$Allele[i] &
                                                    sim.slopes$Slope > true.slopes$Slope[i]))

  }


  multiplot(p1, p2, cols = 1)

  x <- list(true.slopes, sim.slopes)
  names(x) <- c("true.slopes", "simulated.slopes")

  return(x)



}
