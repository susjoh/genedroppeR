#' Summarise genotype frequencies by cohort in a pedigree at a single locus.
#'
#' Summarise genotype frequencies by cohort in a pedigree at a single locus.
#'
#' @param id vector. Individual IDs
#' @param mother vector. Maternal IDs corresponding to id. Missing mothers
#'   should be coded as 0 or NA
#' @param father vector. Paternal IDs corresponding to id. Missing fathers
#'   should be coded as 0 or NA
#' @param cohort vector. Cohort IDs corresponding to id.
#' @param genotype vector. Genotypes IDs corresponding to id.
#' @param na.include boolean. Default FALSE. Indicate if NA cohorts are to be
#'   output in the locus summary.
#' @import plyr
#' @import reshape2
#' @import kinship2
#' @import magrittr
#' @import ggplot2
#' @export

# library(genedroppeR)
#
# data("unicorn")
# head(unicorn)
# table(unicorn$cohort)
# unicorn <- subset(unicorn, cohort %in% 1:10)
#
#
# id       <- unicorn$id
# mother   <- unicorn$mother
# father   <- unicorn$father
# cohort   <- unicorn$cohort
# genotype <- unicorn$genotype
# genotype.delim = ""
#
#   require(plyr)
#   require(reshape2)
#   require(kinship2)
#   require(magrittr)
#   require(ggplot2)
#
#   ped <- data.frame(ID     = id,
#                     MOTHER = mother,
#                     FATHER = father)
#
#   ped$MOTHER[which(ped$MOTHER == 0)] <- NA
#   ped$FATHER[which(ped$FATHER == 0)] <- NA
#
#   if
#   alleles <- as.character(genotype) %>%
#     unique %>%
#     sort %>%
#     strsplit(split = delim) %>%
#     unlist %>%
#     unique
#
#   #~~ Summarise the data by cohort
#
#   suppressMessages(x <- table(ped$cohort, ped$genotype, useNA = "always") %>%
#                      data.frame %>%
#                      dcast(ped$cohort ~ ped$genotype))
#
#   names(x)[1] <- "cohort"
#
#   x$GenoCount <- rowSums(x[,2:(ncol(x)-1)])
#   x$FullCount <- rowSums(x[,2:(ncol(x)-1)])
#   x$PropGenotyped <- x$GenoCount/x$FullCount
#
#   #~~ get allele frequencies by year
#
#   x.freq <- list()
#
#   for(i in alleles){
#
#     atab <- x[,grep(i, genos, value = T)]
#     ahom <- which(names(atab) == paste0(i, delim, i))
#     ahet <- which(names(atab) != paste0(i, delim, i))
#
#     for(j in ahet){
#       atab[,j] <- 0.5*atab[,j]
#     }
#
#     x.freq[[i]] <- rowSums(atab)/x$GenoCount
#
#     rm(atab, ahom, ahet, j)
#
#   }
#
#   x.freq <- data.frame(x.freq)
#
#   names(x.freq) <- paste0("Freq_", names(x.freq))
#
#   x <- cbind(x, x.freq)
#
#   rm(x.freq)
#
#   names(x)[which(names(x) %in% genos)] <- paste0("Geno_", names(x)[which(names(x) %in% genos)])
#
#   x$cohort <- as.numeric(as.character(x$cohort))
#
#   #~~ How many founders?
#
#   if(!is.null(ped)){
#
#     ped$founder <- kindepth(ped[,1], ped[,2], ped[,3])
#     ped$cohort <- cohort
#
#     suppressMessages(x1 <- table(ped$cohort, ped$founder == 0) %>% data.frame %>% dcast(Var1 ~ Var2))
#
#     names(x1) <- c("cohort", "NonFounders", "Founders")
#
#     x1$PropFounders <- x1$Founders/(x1$NonFounders + x1$Founders)
#
#     x1$cohort <- as.character(x1$cohort)
#
#     suppressMessages(x <- join(x, x1))
#   }
#
#   if(na.include == F) if(length(which(is.na(x$cohort))) > 0) x <- subset(x, !is.na(cohort))
#
#   return(x)
#
#
#
# }
