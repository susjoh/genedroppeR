## `genedroppeR`: An R package to conduct single-locus genedrop analyses through pedigrees.

This package conducts "gene-drop" simulations at individual loci through complex pedigree structures to determine if changes in allele frequency are consistent with drift, directional selection, or balancing selection.
### Quick-Start Guide

TO install directly from GitHub:


``` 
library(devtools)
install_github("susjoh/genedroppeR")
library(genedroppeR)
```

The library comes with the `unicorn` dataset from a long-term study on the island of Áiteigin in Scotland (Figure 1). They have been genotyped for three loci: HornSNP, a biallelic locus for horn length; MHC, which characterises multi-allelic variation at the Magic Histocompatibility Locus; and ColourSNP, a SNP responsible for a rare glitter coat polymorphism.

![Figure 1: A typical unicorn.](tutorial/figs/Unicorn.png){width="200"}

Unicorns have a promiscuous mating system, in common with other magical beasts. We can see an example of the pedigree below. The numbers at the side indicate the cohorts of when each unicorn was born. This can be any sequential series of numbers.


