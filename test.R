
library(genedroppeR)
library(dplyr)

data(mhc)


#
id = mhc$id
mother = mhc$mother
father = mhc$father
cohort = mhc$cohort
genotype = mhc$Glitter
nsim = 100
sex = mhc$sex
n_founder_cohorts = 1
fix_founders = T
verbose = T
interval = 10
resample_offspring = F
return_full_results = NULL
remove_founders = TRUE
genotype_delim = ""
multiallelic = F


# library(dplyr)
#
# library(genedroppeR)

data(unicorn)

x <- summary_cohort(id = mhc$id,
               mother = mhc$mother,
               father = mhc$father,
               cohort = mhc$cohort,
               genotype = mhc$Xlinked)

plot_genedrop_cohort(x)


genedrop_obj <- genedrop_snp(id = mhc$id,
                  mother = mhc$mother,
                  father = mhc$father,
                  cohort = mhc$cohort,
                  genotype = mhc$Xlinked,
                  nsim = 100,
                  n_founder_cohorts = 4,
                  fix_founders = T,
                  verbose = T,
                  interval = 1,
                  resample_offspring = F)

summary_genedrop(genedrop_obj)
plot_genedrop(genedrop_obj)

