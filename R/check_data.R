#' Check data before gene-drop analysis.
#'
#' Check data before gene-drop analysis.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id.
#' @param father vector. Paternal IDs corresponding to id.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param genotype vector. Genotypes corresponding to id.
#' @import dplyr
#' @export

check_data <- function(id,
                       mother,
                       father,
                       cohort = NULL,
                       genotype,
                       multiallelic = F){

  #~~ Check that there are no duplicate IDs.

  if (any(as.numeric(names(table(table(id)))) > 1)){
    stop ("Duplicate values in id")
  }

  #~~ If cohort is provided, then check there are no NA's.

  if (!is.null(cohort) & any(is.na(cohort))){
    message("NAs present in cohort - these individuals will be removed.")
  }



  #~~ If genotype is provided, then check the locus is not monomorphic.

  if (!is.null(genotype) & length(table(genotype)) == 1){
    stop ("Locus is monomorphic")
  }

  #~~ If genotype is numeric, then only accept if 0, 1, 2

  if(multiallelic == F){
    if (!is.null(genotype) & is.numeric(genotype)){

      if (any(!na.omit(genotype) %in% 0:2)){
        stop ("Dosage should be coded as 0, 1, 2.")
      }

      if (length(unique(genotype)) == 2){
        message ("This dataset only has two unique genotypes. Please note that this model assumes a value of 1 is the heterozygote")
      }

    } else {

      templev <- as.factor(genotype) %>% levels %>% sort

      #~~ determine which value is heterozygote if templev is shorter than 3

      if (length(templev) == 2){

        templev2 <- sapply(templev, function(x){
          length(unique(strsplit(x, split = "")[[1]]))
        })

        if (templev2[1] == max(templev2)){
          templev <- c("AA", names(templev2))
        }

        if (templev2[1] == templev2[2]){
          templev <- c(names(templev2)[1], "AB", names(templev2)[2])
        }

        rm(templev2)
      }

      genotype <- as.numeric(factor(genotype, levels = templev)) -1
      rm(templev)
    }
  } else {
    genotype <- as.character(genotype)
  }

  # Make a ped object.

  ped <- data.frame(ID     = id,
                    MOTHER = mother,
                    FATHER = father)

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA

  #~~ Add the genotype information.

  if (!is.null(genotype)) ped$genotype <- genotype

  #~~ Add cohort information to ped.

  if (!is.null(cohort)) {

    ped$cohort <- cohort

    #~~ remove anything not cohorted

    if(length(which(is.na(cohort))) > 0){
      message(paste0("Removed ", length(which(is.na(cohort))), " ids with no cohort information."))
    }
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
      stop ("Some parents and offspring are in the same cohort.
            See output for problem lines.")
    }

    rm(x, x1, badids)

  } else {

    ped$cohort <- kindepth(ped$ID, ped$FATHER, ped$MOTHER)
    message(paste0("No cohorts defined. kindepth2::kindepth indicates ", max(ped$cohort), " generations."))

  }

  #~~ If cohort is provided, throw error if any have no genotype information


  x9 <- data.frame(cohort = ped$cohort, genotype = ped$genotype) %>%
    group_by(cohort) %>%
    summarise(n = length(which(!is.na(genotype))))

  if(any(x9$n == 0)){
    x9 <- subset(x9, n == 0)
    stop(paste0("Cohorts ", paste(x9$cohort, collapse = ", "), " have no genotyped individuals."))
  }


  #~~ Get rid of parents that aren't in the IDs.

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA


  return(ped)
}
