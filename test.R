id = unicorn$id
mother = unicorn$mother
father = unicorn$father
cohort = unicorn$cohort
genotype = unicorn$HornSNP
nsim = 100
sex = NULL
n_founder_cohorts = 1
fix_founders = T
verbose = T
interval = 10
resample_offspring = F
return_full_results = NULL

library(dplyr)

library(genedroppeR)

data(unicorn)

x <- summary_cohort(id = unicorn$id,
               mother = unicorn$mother,
               father = unicorn$father,
               cohort = unicorn$cohort,
               genotype = unicorn$HornSNP)
summary(x)
plot(x)

genedrop_obj <- genedrop_snp(id = unicorn$id,
                  mother = unicorn$mother,
                  father = unicorn$father,
                  cohort = unicorn$cohort,
                  genotype = unicorn$HornSNP,
                  nsim = 10,
                  n_founder_cohorts = 4,
                  fix_founders = T,
                  verbose = T,
                  interval = 1,
                  resample_offspring = F)

summary(x)
plot(x)

