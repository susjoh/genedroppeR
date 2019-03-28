
library(genedroppeR)

data("unicorn")

x <- locus_summary(id  = unicorn$id,
              mother   = unicorn$mother,
              father   = unicorn$father,
              cohort   = unicorn$cohort,
              genotype = unicorn$genotype)
x

plot_locus_summary(x)

x

sims <- genedrop(id       = unicorn$id,
                 mother   = unicorn$mother,
                 father   = unicorn$father,
                 cohort   = unicorn$cohort,
                 genotype = unicorn$genotype,
                 nsim = 100,
                 n_founder_cohorts = 7,
                 fix_founders = T)

plot_genedrop_slopes(sims)
str(sims)

