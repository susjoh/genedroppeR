#' Plot locus summary output
#'
#' Plot locus summary output to investigate proportion of IDs genotyped,
#' proportion of founder individuals and the total ID count per cohort.
#'
#' @param x Output of locus_summary() function
#' @import ggplot2
#' @export
#'

plot_locus_summary <- function(x){

  #~~ What proportion of IDs are genotyped from year to year?

  p1 <- ggplot(x, aes(cohort, PropGenotyped)) +
    geom_line(colour = "red") +
    labs(x = "Cohort", y = "Proportion of IDs genotyped")

  #~~ How many genotyped individuals and how many individuals in total?

  x.temp <- melt(x[,c("cohort", "GenoCount", "FullCount")], id.vars = "cohort")

  p2 <- ggplot(x.temp, aes(cohort, value, colour = variable)) +
    geom_line() +
    labs(x = "Cohort", y = "Total ID Count per Cohort") +
    theme(legend.position = "top")

  #~~ What proportion of IDs are pedigree founders from year to year?

  p3 <- ggplot(x, aes(cohort, PropFounders)) +
    geom_line(colour = "red") +
    labs(x = "Cohort", y = "Proportion of IDs that are founders")

  print(multiplot(p1, p2, p3, cols = 2))



}
