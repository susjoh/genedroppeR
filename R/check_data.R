#' `check_data()`: Check data formatting before genedrop analysis.
#'
#' This function is used within `genedrop_...` functions to check and format the
#' data for downstream genedrop analysis. Returns an object `ped` that is
#' correctly formatted for genedrop analysis. Note - this function is
#' automatically run in all `genedrop_...` functions.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id. Missing values should
#'   be 0 or NA.
#' @param father vector. Paternal IDs corresponding to id. Missing values should
#'   be 0 or NA.
#' @param cohort vector (optional). Cohort number (e.g. birth year)
#'   corresponding to the id.
#' @param genotype vector. Genotypes corresponding to id.
#' @param sex vector (optional). Sexes corresponding to id. 1 is the homogametic
#'   sex (e.g. XY, ZW) and 2 is the heterogametic sex (e.g. XX, ZZ). Specifying
#'   `sex` will force strict correspondence to maternal and paternal links. This
#'   vector only needs to be specified for simulations on sex-linked SNPs.
#' @param multiallelic boolean. Default = FALSE
#' @import dplyr
#' @export

check_data <- function(id,
                       mother,
                       father,
                       cohort = NULL,
                       genotype,
                       sex = NULL,
                       multiallelic = FALSE){

  ID = MOTHER = FATHER = NULL

  #

  # Check that there are no duplicate IDs.

  if (any(as.numeric(names(table(table(id)))) > 1)){
    stop ("Duplicate values in id")
  }

  # If cohort is provided, then check there are no NA's and that it is an integer

  if (!is.null(cohort) & !is.integer(cohort)){
    stop("Cohort must be an integer")
  }

  if (!is.null(cohort) & any(is.na(cohort))){
    message("NAs present in cohort - these individuals will be removed.")
  }

  # If genotype is provided, then check the locus is not monomorphic.

  if (!is.null(genotype) & length(table(genotype)) == 1){
    stop ("Locus is monomorphic")
  }

  # If genotype is numeric, then only accept if 0, 1, 2

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

      # determine which value is heterozygote if templev is shorter than 3

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

  ped <- data.frame(ID     = as.character(id),
                    MOTHER = as.character(mother),
                    FATHER = as.character(father))

  ped$MOTHER[which(ped$MOTHER == 0)] <- NA
  ped$FATHER[which(ped$FATHER == 0)] <- NA

  # Add the genotype information.

  if (!is.null(genotype)) ped$genotype <- genotype

  # Add the sex information.

  if (!is.null(sex)) ped$sex <- sex

  # Add cohort information to ped.

  if (!is.null(cohort)) {

    ped$cohort <- cohort

    # remove anything not cohorted

    if(length(which(is.na(cohort))) > 0){
      message(paste0("Removed ", length(which(is.na(cohort))), " ids with no cohort information."))
    }
    ped <- filter(ped, !is.na(cohort))

    # Check that parents and offspring are not in the same cohort.

    x <- select(ped, ID, MOTHER, FATHER, cohort)

    x1 <- select(ped, ID, cohort)

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
    message(paste0("No cohorts defined - cohorts determined based on pedigree depth. There are ", max(ped$cohort)+1, " generations."))

  }

  # If cohort is provided, throw error if any have no genotype information

  x9 <- data.frame(cohort = ped$cohort, genotype = ped$genotype) %>%
    group_by(cohort) %>%
    summarise(n = length(which(!is.na(genotype))))

  if(any(x9$n == 0)){
    x9 <- select(x9, n == 0)
    stop(paste0("Cohorts ", paste(x9$cohort, collapse = ", "), " have no genotyped individuals."))
  }

  # Blank out parents that aren't in the IDs (no genotype information).

  ped$MOTHER[which(!is.na(ped$MOTHER) & !ped$MOTHER %in% ped$ID)] <- NA
  ped$FATHER[which(!is.na(ped$FATHER) & !ped$FATHER %in% ped$ID)] <- NA

  # Check all the sexes are okay re: heterogametic sex

  if(!is.null(sex)){

    if(any(is.na(ped$sex))){
      stop("NA values in sex - sex must be specified for sex-linked loci.")
    }

    if(length(unique(ped$sex[which(ped$ID %in% ped$MOTHER)])) > 1){
      stop("Some mothers have been assigned as males - please check sex designations.")
    }

    if(length(unique(ped$sex[which(ped$ID %in% ped$FATHER)])) > 1){
      stop("Some fathers have been assigned as females - please check sex designations.")
    }

    if(any(ped$ID[which(ped$sex == 1)] %in% ped$MOTHER) || any(ped$ID[which(ped$sex == 2)] %in% ped$FATHER)){
      message("Sex-linked models will assume a ZW system")
    } else {
      message("Sex-linked models will assume a XY system")
    }

  }

  ped
}
