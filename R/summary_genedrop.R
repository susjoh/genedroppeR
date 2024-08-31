#' Generic summary function for class "genedroppeR"
#'
#' @param genedrop_obj An object of class "genedroppeR" as output by functions `genedrop_...()`.
#' @param ... Further arguments to be passed
#' @exportS3Method base::summary

summary.genedroppeR <- function(genedrop_obj, ...){
  cat("\n")
  cat("*** Summary: genedroppeR analysis ***")
  cat("\n\n")
  cat(paste0("Model parameters: cohorts = ", length(unique(genedrop_obj$observed_frequencies$Cohort)),
             ", Simulations = ", max(genedrop_obj$simulated_frequencies$Simulation),
             ", Founder cohorts = ", genedrop_obj$n_founder_cohorts, "."))
  cat("\n")
  cat(paste0("Founder genotypes ", ifelse(genedrop_obj$fix_founders, "fixed", "not fixed"),
             ", Offspring ", ifelse(genedrop_obj$resample_offspring, "resampled", "not resampled"), "."))
  cat("\n")
  cat(paste0("Founder cohorts ", ifelse(genedrop_obj$remove_founders, "have been removed from the below analyses (Default).",
                                        "have not been removed from the below analyses.")))
  cat("\n\n")
  data.frame(genedrop_obj$results)

}

#' Generic summary function for class "genedroppeR_cohort"
#' @param cohort_obj An object of class "cohort_obj" as output by
#'   function `summary_cohort()`.
#' @param ... Further arguments to be passed
#' @exportS3Method base::summary

summary.genedroppeR_cohort <- function(cohort_obj, ...){
  class(cohort_obj) <- NULL
  data.frame(cohort_obj)
}


#' Summarise allele frequencies from a genedrop object.
#'
#' Summarise observed and simulated allele frequencies by simulation and cohort
#' for a genedrop object. This is only called within specific functions and is
#' not for general use.
#'
#' @param genedrop_obj An object of class "genedroppeR" as output by functions
#'   `genedrop_...()`.
#' @param genotype_delim char. Default = ''. Delimiter character for genotypes.
#' @export


summarise_genedrop <- function(genedrop_obj, genotype_delim = ''){

  True.Geno = Simulation = cohort = NULL

  genedrop_obj <- filter(genedrop_obj, !is.na(True.Geno))

  if(is.numeric(genedrop_obj$Simulated.Geno)){

    x1 <- melt(tapply(genedrop_obj$Simulated.Geno, list(genedrop_obj$cohort, genedrop_obj$Simulation), sum))
    x2 <- melt(tapply(genedrop_obj$Simulated.Geno, list(genedrop_obj$cohort, genedrop_obj$Simulation), length))
    names(x2)[3] <- "Count"

    suppressMessages(x1 <- left_join(x1, x2))
    names(x1) <- c("Cohort", "Simulation", "Sum", "Count")
    x1$p <- 1 - (x1$Sum/(2*x1$Count))

    x1$Sum <- NULL

    x3 <- filter(genedrop_obj, Simulation == 1) %>% select(True.Geno, cohort) %>% na.omit
    x3 <- data.frame(Sum = tapply(x3$True.Geno, x3$cohort, sum),
                     Count = tapply(x3$True.Geno, x3$cohort, length))

    x3$p <- 1 - (x3$Sum/(2*x3$Count))
    x3$Simulation <- 0
    x3$Cohort <- as.numeric(row.names(x3))
    x3$Sum <- NULL

    x3 <- x3[,c("Cohort", "Simulation", "Count", "p")]

  } else {

    # simulated frequencies

    x1 <- melt(genedrop_obj[,c("cohort", "Simulation", "Mum.Allele", "Dad.Allele")], id.vars = c("cohort", "Simulation"))
    x2a <- data.frame(table(x1$cohort, x1$Simulation, x1$value))
    x2b <- data.frame(table(x1$cohort, x1$Simulation))
    names(x2b)[3] <- "Count"
    suppressMessages(x1 <- left_join(x2a, x2b))

    x1$Freq <- x1$Freq/x1$Count

    names(x1) <- c("Cohort", "Simulation", "Allele", "p", "Count")

    x1 <- x1[,c("Cohort", "Simulation", "Count", "p", "Allele")]

    rm(x2a, x2b)

    # observed frequencies

    x3 <- filter(genedrop_obj, Simulation == 1) %>% select(True.Geno, cohort) %>% na.omit

    x3$Allele1 <- sapply(x3$True.Geno, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][1])
    x3$Allele2 <- sapply(x3$True.Geno, function(foo) strsplit(foo, split = genotype_delim, fixed = T)[[1]][2])

    x4 <- melt(x3[,c("cohort", "Allele1", "Allele2")], id.vars = c("cohort"))

    x2a <- data.frame(table(x4$cohort, x4$value))
    x2b <- data.frame(table(x4$cohort))
    names(x2b)[2] <- "Count"
    suppressMessages(x3 <- left_join(x2a, x2b))

    x3$Freq <- x3$Freq/x3$Count
    x3$Simulation <- 0

    names(x3) <- c("Cohort", "Allele", "p", "Count", "Simulation")

    x3 <- x3[,c("Cohort", "Simulation", "Count", "p", "Allele")]

    rm(x2a, x2b, x4)

  }

  list(observed_frequencies = x3,
       simulated_frequencies = x1)

}
