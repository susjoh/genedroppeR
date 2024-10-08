#' `resample_offspring_func()`: Reassign offspring across parents in each
#' cohort.
#'
#' This function will reassign offspring to each parent to ameliorate an
#' occasional criticism of the genedrop method, that the pedigree is a product
#' of fitness (i.e. "fitter" parents may have "fitter" offspring). This function
#' reassigns offspring to each breeding parent within a cohort and is  run in
#' `genedrop_...` functions when specifying `resample_offspring = T`. This will
#' preserve the distribution of offspring in mothers and fathers. At present,
#' non-breeding individuals are not sampled as parents as it cannot be
#' determined if they are available to mate based on pedigree information alone.
#'
#' BE AWARE: using this function means that a direct comparison to the true
#' pedigree is no longer possible - it is up to the user to make a philosophical
#' choice on what to do. More discussion is provided on this in the package
#' vignette.
#'
#' @param ped A pedigree object as generated by check_data
#' @returns pedigree with resampled offspring
#' @noRd


resample_offspring_func <- function(ped) {
  cohort <- MOTHER <- FATHER <- NULL

  lu_mum <- NULL
  lu_dad <- NULL

  for (h in 1:max(ped$cohort)) {
    coparents <- filter(ped, cohort == h)$ID

    lookup_mum <- data.frame(MOTHER = sort(unique(filter(ped, MOTHER %in% coparents)$MOTHER)))

    if (nrow(lookup_mum) > 0) {
      if (nrow(lookup_mum) > 1) {
        lookup_mum$new.MOTHER <- sample(lookup_mum$MOTHER, replace = F)
      } else {
        lookup_mum$new.MOTHER <- lookup_mum$MOTHER
      }

      lu_mum <- rbind(lu_mum, lookup_mum)
    }

    lookup_dad <- data.frame(FATHER = sort(unique(filter(ped, FATHER %in% coparents)$FATHER)))

    if (nrow(lookup_dad) > 0) {
      if (nrow(lookup_dad) > 1) {
        lookup_dad$new.FATHER <- sample(lookup_dad$FATHER, replace = F)
      } else {
        lookup_dad$FATHER <- lookup_dad$FATHER
      }

      lu_dad <- rbind(lu_dad, lookup_dad)
    }

    rm(lookup_mum, lookup_dad, coparents)
  }

  suppressMessages(ped <- left_join(ped, lu_mum))
  suppressMessages(ped <- left_join(ped, lu_dad))

  ped$MOTHER <- ped$new.MOTHER
  ped$FATHER <- ped$new.FATHER

  ped
}
