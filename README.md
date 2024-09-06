[![pkgcheck](https://github.com/susjoh/genedroppeR/workflows/pkgcheck/badge.svg)](https://github.com/susjoh/genedroppeR/actions?query=workflow%3Apkgcheck) [![R-CMD-check](https://github.com/susjoh/genedroppeR/workflows/R-CMD-check/badge.svg)](https://github.com/susjoh/genedroppeR/actions?query=workflow%3AR-CMD-check)\

## `genedroppeR`: An R package to conduct single-locus genedrop analyses through pedigrees.

This package conducts "genedrop" simulations at individual loci through complex pedigree structures to determine if changes in allele frequency are consistent with e.g. drift, directional selection, or balancing selection.

### Quick-Start Guide

To install directly from GitHub:

```         
library(remotes)
install_github("susjoh/genedroppeR")
```

### Case study: the unicorns of Áiteigin

![Figure 1: A typical unicorn.](tutorial/figs/Unicorn.png)

The library comes with the `unicorn` dataset from a long-term study on the island of Áiteigin in Scotland (Figure 1). This can be explored by loading `data(unicorn)` They have been genotyped for four loci:

1.  Horns, a biallelic locus for horn length;
2.  MHC, which characterises multi-allelic variation at the Magic Histocompatibility Locus;
3.  Glitter, a biallelic locus responsible for a rare glitter coat polymorphism;
4.  Xlinked, a generic X-linked SNP.

```         
> head(unicorn)
    id mother father cohort sex Horns MHC Xlinked Glitter
1 -136      0      0   2001   2    GG  BD      GA      GG
2 -135      0      0   2001   2    AA  GG      AA      AG
3 -131      0   4552   2001   2    AA  FH      GG      GG
4 -130      0   4554   2001   2    GA  EH      GA      GG
5 -129      0   4555   2001   2    GA  DG      GA      GG
6 -128      0   4556   2001   2    GA  BH      GG      GG
```

Unicorns have a polygamous mating system, in common with other magical beasts. We can see an example of the pedigree below. The numbers at the side indicate the cohorts of when each unicorn was born:

![Figure 2: A unicorn pedigree.](tutorial/figs/pedigree.png)

### Step 1: Summarise the observed data.

We can summarise this data in terms of the allele frequency changes using the function `summary_cohort()`. Let's start with the HornSNP locus.

```         
unicorn_summ <- summary_cohort(id = unicorn$id,
                               mother = unicorn$mother,
                               father = unicorn$father,
                               cohort = unicorn$cohort,
                               genotype = unicorn$Horns)

plot_genedrop_cohort(unicorn_summ)
```

This output provides the allele frequencies for the A and G alleles per cohort, as well as counts of (genotyped) individuals per cohort, the proportion of genotyped individuals, founders and non-founders, and the proportion of founders per cohort. A founder is defined as an individual where neither the mother or father is known. The frequency of the A allele, which confers a larger horn, seems to be increasing in the population:

![Figure 3: Output from `summary_cohort()` and `plot_genedrop_cohort()`.](tutorial/figs/summary_cohort.png)

However, given the complex structure of the pedigree, it may be that this increase in allele frequency is more likely to be due to drift than selection. One approach to model this is to simulate allele frequency changes given a pedigree of the same structure, and see if the observed allele frequency change is significantly different from what we would expect by chance. This is the motivation for conducting a "genedrop" analysis, where we can model this based on the observed pedigree structure.

#### Example 1: The Horns locus

In early cohorts, many individuals are likely to be founders. To run the genedrop, we select some cohorts at the start to be the "founder" cohorts, where allele frequencies are sampled from the observed frequencies in each year. In Figure 3, the number of cohorts declines rapidly until around 5 cohorts (1989 - 1993), so lets define the first 4 cohorts as founder cohorts, and the rest as simulated cohorts i.e. we will investigate how allele frequencies change in the second 26 year period from 1993 to 2018.

As *Horns* is a biallelic locus, we use the function `genedrop_snp()`. We define the `id`, `mother`, `father`, `cohort` and `genotype` using columns from the `unicorn` dataset. We then define `nsim`, the number of simulations to run and the number of founder cohorts. We can also select `fix_founders` which keeps the founder genotypes fixed in genotyped individuals, i.e. how allele frequencies would change from an almost identical starting point.

```         
genedrop_obj <- genedrop_snp(
   id = unicorn$id,
   mother = unicorn$mother,
   father = unicorn$father,
   cohort = unicorn$cohort,
   genotype = unicorn$Horns,
   nsim = 1000,
   n_founder_cohorts = 4,
   fix_founders = TRUE)
```

We can then get the summary results by running:

```         
summary_genedrop(genedrop_obj)

> summary_genedrop(genedrop_obj)

*** Summary: genedroppeR analysis ***

Model parameters: cohorts = 30, Simulations = 1000, Founder cohorts = 4.
Founder genotypes fixed, Offspring not resampled.
Founder cohorts have been removed from the below analyses (Default).

            Analysis Allele    Estimate Estimates.Lower Estimates.Higher Simulations
1  Cumulative Change      p 0.709569448             358              642        1000
2 Directional Change      p 0.002192312             869              131        1000

plot_genedrop(genedrop_obj)
```

![Figure 4: Output from `plot_genedrop()`.](tutorial/figs/genedrop_horns.png)
