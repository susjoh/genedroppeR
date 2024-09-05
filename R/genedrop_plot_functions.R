#' `plot_genedrop()`: Plot a summary graph of `genedrop_...()` outputs.
#'
#' @param genedrop_obj Gene-Drop object
#' @param analysis_to_plot Default = "all". Options are "frequency",
#'   "directional" and "cumulative". This will output graphs of allele frequency
#'   changes, directional change, and cumulative change. If the locus is
#'   multi-allelic and the user specifies the default option ("all"), the user
#'   will be prompted to press <Enter> to advance through each option.
#' @param sim_alpha Default = 0.3. Alpha (transparency) value for plotted lines.
#' @param obs_line_col Default = "red". Line colour to use for the observed
#'   data.
#' @examples
#' data(unicorn)
#' genedrop_obj <- genedrop_snp(
#'   id = unicorn$id,
#'   mother = unicorn$mother,
#'   father = unicorn$father,
#'   cohort = unicorn$cohort,
#'   genotype = unicorn$Horns,
#'   nsim = 100,
#'   n_founder_cohorts = 4,
#'   fix_founders = TRUE,
#'   verbose = TRUE,
#'   interval = 10,
#'   resample_offspring = FALSE
#' )
#' plot_genedrop(genedrop_obj)
#'
#' @export
#'

plot_genedrop <- function(genedrop_obj,
                          analysis_to_plot = "all",
                          sim_alpha = 0.3,
                          obs_line_col = "red") {
  Cohort <- Simulation <- p <- Estimate <- NULL


  alleles <- sort(unique(genedrop_obj$simulated_frequencies$Allele))

  p1 <- ggplot(genedrop_obj$simulated_frequencies, aes(Cohort, p, group = Simulation)) +
    geom_line(alpha = sim_alpha) +
    geom_line(data = genedrop_obj$observed_frequencies, aes(Cohort, p), col = obs_line_col) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    geom_vline(xintercept = min(genedrop_obj$observed_frequencies$Cohort) + genedrop_obj$n_founder_cohorts - 0.5, linetype = "dashed") +
    ggtitle(paste0("Allele Frequency Changes: Nsim = ", max(genedrop_obj$simulated_frequencies$Simulation))) +
    facet_wrap(~Allele)

  p2 <- ggplot(genedrop_obj$slopes$sim.slopes, aes(Estimate)) +
    geom_histogram(col = "grey") +
    facet_wrap(~Allele) +
    geom_vline(data = genedrop_obj$slopes$true.slopes, aes(xintercept = Estimate), col = obs_line_col) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    ggtitle(paste0("Distribution of Regression Slopes: Nsim = ", max(genedrop_obj$simulated_frequencies$Simulation)))

  p3 <- ggplot(genedrop_obj$cumulative_change$sim.changes, aes(Estimate)) +
    geom_histogram(col = "grey") +
    facet_wrap(~Allele) +
    geom_vline(data = genedrop_obj$cumulative_change$true.changes, aes(xintercept = Estimate), col = obs_line_col) +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    ggtitle(paste0("Distribution of Cumulative Change: Nsim = ", max(genedrop_obj$simulated_frequencies$Simulation)))

  if (analysis_to_plot == "all") {
    if (length(alleles) == 1) {
      suppressMessages(gridExtra::grid.arrange(p1, p2, p3, layout_matrix = rbind(
        c(1, 1),
        c(2, 3)
      )))
    } else {
      print(p1)
      message("Press <Enter> to continue...")
      readline() # Waits for the user to press enter
      print(p2)
      message("Press <Enter> to continue...")
      readline() # Waits for the user to press enter
      print(p3)
    }
  }

  if (analysis_to_plot == "frequency") print(p1)
  if (analysis_to_plot == "directional") print(p2)
  if (analysis_to_plot == "cumulative") print(p3)
}

#' `plot_genedrop_cohort()`: Plot a summary graph of an object produced by
#' `summary_cohort()`. This will plot the observed temporal dynamics of the
#' locis, as well as the proprtion of founder individuals and genotyped
#' individuals within each cohort.
#'
#' @param cohort_obj data.frame created by `summary_cohort()`
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
#'


plot_genedrop_cohort <- function(cohort_obj) {
  cohort <- PropFounders <- PropGenotyped <- variable <- value <- NULL

  x <- cohort_obj
  x <- subset(x, cohort != "missing")
  x$cohort <- as.numeric(as.character(x$cohort))

  # Proportion founders and genotyped

  p2 <- ggplot(x, aes(cohort, PropFounders)) +
    geom_line() +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    ggtitle("Proportion of Founders")

  p3 <- ggplot(x, aes(cohort, PropGenotyped)) +
    geom_line() +
    theme_bw() +
    theme(plot.title = element_text(size = 10)) +
    ggtitle("Proportion of Cohort Genotyped")

  # Temporal trend


  x1 <- melt(x, id.vars = c(
    "cohort", "GenoCount", "FullCount", "PropGenotyped",
    "NonFounders", "Founders", "PropFounders"
  ))

  x1$variable <- as.character(x1$variable)

  alleles <- unique(sort(x1$variable))


  if (length(alleles) < 3) {
    x1 <- subset(x1, variable == alleles[1])

    p1 <- ggplot(x1, aes(cohort, value)) +
      geom_line() +
      stat_smooth(method = "lm") +
      facet_wrap(~variable) +
      theme_bw() +
      theme(plot.title = element_text(size = 10)) +
      ggtitle("Temporal dynamics of alleles")
  } else {
    p1 <- ggplot(x1, aes(cohort, value)) +
      geom_line() +
      stat_smooth(method = "lm") +
      facet_wrap(~variable) +
      theme_bw() +
      theme(plot.title = element_text(size = 10)) +
      ggtitle("Temporal dynamics of alleles")
  }

  suppressMessages({
    if (length(alleles) < 3) {
      gridExtra::grid.arrange(p1, p2, p3, layout_matrix = rbind(
        c(1, 1),
        c(2, 3)
      ))
    } else {
      print(p1)
      message("Press <Enter> to continue...")
      readline() # Waits for the user to press enter
      gridExtra::grid.arrange(p2, p3, layout_matrix = rbind(
        c(1),
        c(2)
      ))
    }
  })
}
