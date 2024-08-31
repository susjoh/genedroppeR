
library(genedroppeR)
library(dplyr)

data(unicorn)


#
id = unicorn$id
mother = unicorn$mother
father = unicorn$father
cohort = unicorn$cohort
genotype = unicorn$MHC
genotype_delim = ""
nsim = 100
n_founder_cohorts = 1
fix_founders = T
verbose = T
interval = 10
resample_offspring = F
return_full_results = NULL
remove_founders = TRUE
multiallelic = T


# library(dplyr)
#
# library(genedroppeR)

data(unicorn)

unicorn <- subset(unicorn, cohort < 2010)

x <- summary_cohort(id = unicorn$id,
               mother = unicorn$mother,
               father = unicorn$father,
               cohort = unicorn$cohort,
               genotype = unicorn$MHC)

plot_genedrop_cohort(x)


genedrop_obj <- genedrop_multi(id = unicorn$id,
                  mother = unicorn$mother,
                  father = unicorn$father,
                  cohort = unicorn$cohort,
                  genotype = unicorn$MHC,
                  nsim = 100,
                  n_founder_cohorts = 4,
                  fix_founders = T,
                  verbose = T,
                  interval = 1,
                  resample_offspring = F)

summary_genedrop(genedrop_obj)
plot_genedrop(genedrop_obj)

