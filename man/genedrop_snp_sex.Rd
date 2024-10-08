% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genedrop_snp_sexlinked.R
\name{genedrop_snp_sex}
\alias{genedrop_snp_sex}
\title{`genedrop_snp_sexlinked()`: Conduct a genedrop simulation for a sex-linked
biallelic locus.}
\usage{
genedrop_snp_sex(
  id,
  mother,
  father,
  cohort = NULL,
  genotype,
  sex,
  nsim,
  n_founder_cohorts = 1,
  fix_founders = TRUE,
  verbose = TRUE,
  interval = 100,
  resample_offspring = FALSE,
  remove_founders = TRUE,
  return_full_results = NULL
)
}
\arguments{
\item{id}{vector. Individual IDs}

\item{mother}{vector. Maternal IDs corresponding to id.}

\item{father}{vector. Paternal IDs corresponding to id.}

\item{cohort}{vector (optional). Cohort number (e.g. birth year)
corresponding to the id.}

\item{genotype}{vector. Genotypes corresponding to id.}

\item{sex}{vector. Sexes corresponding to id. 1 is the homogametic sex (e.g.
XY, ZW) and 2 is the heterogametic sex (e.g. XX, ZZ)}

\item{nsim}{integer. Number of genedrop simulations to run.}

\item{n_founder_cohorts}{integer. The number of cohorts at the top of the
pedigree that will sample from the true allele frequencies (these are
defined as "sampled"). All cohorts following these ones are "simulated" and
are used for comparisons of changes in allele frequency.}

\item{fix_founders}{logical. Default = TRUE. Determines whether individuals
in founder cohorts should be given their true recorded genotypes (if
known). For individuals with no known genotype, their genotypes are sampled
based on the observed cohort allele frequency. If FALSE, then all IDs are
sampled based on the cohort allele frequencies.}

\item{verbose}{logical. Default = TRUE. Output the progress of the run.}

\item{interval}{integer. Default 100. Output progress every 100 simulations.}

\item{resample_offspring}{logical. Default = FALSE. If FALSE, the same
pedigree structure as the observed pedigree is used. If TRUE, then
offspring are resampled across parents in each cohort. This is to remove
any potential signal where prolific individuals tend to have prolific
offspring, but will also mean that pedigrees are not directly comparable.}

\item{remove_founders}{Default = TRUE. If TRUE, then the founder cohorts will
be removed from calculations of directional and cumulative change.}

\item{return_full_results}{Default = NULL. This will also output tables of
all individually simulated genotypes.}
}
\value{
an output object of class "genedroppeR"
}
\description{
This function conducts a genedrop simulation for a single, sex-linked
bi-allelic locus (e.g. a, X- or Z-linked SNP). For autosomal loci, use
`genedrop_snp()`. Before running this function, users should first summarise
and visualise their data using `summary_cohort()` to determine an appropriate
value for `n_founder_cohorts`. This function will return an object that
contains the cohort allele frequences in the observed and simulated datasets.
Overall results of directional and balancing selection can be observed using
`summary()`. For more detail on specifying model parameters, please consult
the tutorial at https://github.com/susjoh/genedroppeR.
}
\examples{
data(unicorn)
sub_unicorn <- subset(unicorn, sex \%in\% c(1, 2))
genedrop_obj <- genedrop_snp_sex(
  id = sub_unicorn$id,
  mother = sub_unicorn$mother,
  father = sub_unicorn$father,
  cohort = sub_unicorn$cohort,
  genotype = sub_unicorn$Xlinked,
  sex = sub_unicorn$sex,
  nsim = 100,
  n_founder_cohorts = 4,
  fix_founders = TRUE,
  verbose = TRUE,
  interval = 10
)
summary_genedrop(genedrop_obj)
plot_genedrop(genedrop_obj)
}
