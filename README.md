[![pkgcheck](https://github.com/susjoh/genedroppeR/workflows/pkgcheck/badge.svg)](https://github.com/susjoh/genedroppeR/actions?query=workflow%3Apkgcheck)
[![R-CMD-check](https://github.com/susjoh/genedroppeR/workflows/R-CMD-check/badge.svg)](https://github.com/susjoh/genedroppeR/actions?query=workflow%3AR-CMD-check)
## `genedroppeR`: An R package to conduct single-locus genedrop analyses through pedigrees.

This package conducts "genedrop" simulations at individual loci through complex pedigree structures to determine if changes in allele frequency are consistent with e.g. drift, directional selection, or balancing selection.

### Quick-Start Guide

To install directly from GitHub:

```         
library(remotes)
install_github("susjoh/genedroppeR")
```

### Case study: the unicorns of Áiteigin

The library comes with the `unicorn` dataset from a long-term study on the island of Áiteigin in Scotland (Figure 1). This can be explored by loading `data(unicorn)` They have been genotyped for four loci:

1. Horns, a biallelic locus for horn length;
2. MHC, which characterises multi-allelic variation at the Magic Histocompatibility Locus;
3. Glitter, a biallelic locus responsible for a rare glitter coat polymorphism;
4. Xlinked, a generic X-linked SNP.

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

![Figure 1: A typical unicorn.](tutorial/figs/Unicorn.png)

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

plot_cohort_summary(unicorn_summ)

```

This output provides the allele frequencies for the A and G alleles per cohort, as well as a count of the number of genotyped individuals and the total number of individuals per cohort, the proportion of genotyped individuals, the numbers of non founders and founders, and the proportion of founders per cohort. A founder is defined as an individual where neither the mother or father is known.

The frequency of the A allele, which confers a larger horn, seems to be increasing in the population:


#### Example 1: The Horns locus

