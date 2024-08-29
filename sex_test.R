library(dplyr)
library(genedroppeR)

data(unicorn)

haplo.frame$SexGeno <- apply(haplo.frame[,c("Hom.Parent.Allele", "Het.Parent.Allele")], 1, sum, na.rm = T)
haplo.frame$SexGeno[which(haplo.frame$sex == 1 & haplo.frame$SexGeno == 1)] <- 2
haplo.frame$SexGeno <- ifelse(haplo.frame$SexGeno == 0, "AA",
                              ifelse(haplo.frame$SexGeno == 1, "AB", "BB"))

unicorn$SexLinkedGeno <- haplo.frame$SexGeno

unicorn$SexLinkedGenoWithError <- unicorn$SexLinkedGeno

unicorn$SexLinkedGenoWithError[sample(1:nrow(unicorn), size = 446)] <- NA

unicorn$sex <- haplo.frame$sex

head(unicorn)
unicorn <- unicorn[,c(1, 3, 2, 4, 10, 5:9)]


write.table(unicorn, "data/unicorn.txt", row.names = F, sep = "\t", quote = F)

setwd("test_files/")

unicorn$MHC <- NULL
unicorn$Pheno <- -9
 x <- function(foo) paste(strsplit(foo, split = "")[[1]], collapse = "\t")
unicorn$HornSNP <- sapply(as.character(unicorn$HornSNP), x)
unicorn$ColourSNP  <- sapply(as.character(unicorn$ColourSNP ), x)
unicorn$SexLinkedGeno  <- sapply(as.character(unicorn$SexLinkedGeno ), x)
unicorn$SexLinkedGenoWithError  <- sapply(as.character(unicorn$SexLinkedGenoWithError ), x)
unicorn <- unicorn[,c(1:5, 10, 6:9)]
unicorn[is.na(unicorn)] <- 0

unicorn$HornSNP[which(unicorn$HornSNP == "NA")] <- "0\t0"
unicorn$ColourSNP[which(unicorn$ColourSNP == "NA")] <- "0\t0"
unicorn$SexLinkedGeno[which(unicorn$SexLinkedGeno == "NA")] <- "0\t0"
unicorn$SexLinkedGenoWithError[which(unicorn$SexLinkedGenoWithError == "NA")] <- "0\t0"
unicorn$cohort <- NULL

unicorn <- cbind(Fid = 1, unicorn)

write.table(unicorn, "unicorn.ped", row.names = F, col.names = F, quote = F)

system("plink --file unicorn --mendel")
