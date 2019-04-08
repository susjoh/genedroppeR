
genedrop_summary <- function(genedrop_object){
  
  library(reshape2)
  
  genedrop_object <- subset(genedrop_object, !is.na(True.Geno))
  x1 <- melt(tapply(genedrop_object$Simulated.Geno, list(genedrop_object$cohort, genedrop_object$Simulation), sum))
  x2 <- melt(tapply(genedrop_object$Simulated.Geno, list(genedrop_object$cohort, genedrop_object$Simulation), length))
  names(x2)[3] <- "Count"
  
  x1 <- join(x1, x2)
  head(x1)
  names(x1) <- c("Cohort", "Simulation", "Sum", "Count")
  x1$p <- 1 - (x1$Sum/(2*x1$Count))
  
  x1$Sum <- NULL
  
  x3 <- subset(x, Simulation == 1) %>% subset(select = c(True.Geno, cohort)) %>% na.omit
  x3 <- data.frame(Sum = tapply(x3$True.Geno, x3$cohort, sum),
                   Count = tapply(x3$True.Geno, x3$cohort, length))
  
  x3$p <- 1 - (x3$Sum/(2*x3$Count))
  x3$Simulation <- 0
  x3$Cohort <- as.numeric(row.names(x3))
  x3$Sum <- NULL
  
  head(x1)
  head(x3)
  
  x3 <- x3[,c("Cohort", "Simulation", "Count", "p")]
  
  return(list("observed_frequencies" = x3,
              "simulated_frequencies" = x1))
  
}