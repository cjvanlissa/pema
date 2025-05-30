% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate_smd.R
\name{simulate_smd}
\alias{simulate_smd}
\title{Simulates a meta-analytic dataset}
\usage{
simulate_smd(
  k_train = 20,
  k_test = 100,
  mean_n = 40,
  es = 0.5,
  tau2 = 0.04,
  alpha = 0,
  moderators = 5,
  distribution = "normal",
  model = "es * x[, 1]"
)
}
\arguments{
\item{k_train}{Atomic integer. The number of studies in the training dataset.
Defaults to 20.}

\item{k_test}{Atomic integer. The number of studies in the testing dataset.
Defaults to 100.}

\item{mean_n}{Atomic integer. The mean sample size of each simulated study in
the meta-analytic dataset. Defaults to \code{40}. For each simulated study, the
sample size n is randomly drawn from a normal distribution with mean mean_n,
and sd mean_n/3.}

\item{es}{Atomic numeric vector. The effect size, also known as beta, used in
the model statement. Defaults to \code{.5}.}

\item{tau2}{Atomic numeric vector. The residual heterogeneity. For a range of
realistic values encountered in psychological research, see Van Erp,
Verhagen, Grasman, & Wagenmakers, 2017. Defaults to \code{0.04}.}

\item{alpha}{Vector of slant parameters, passed to \link[sn:dsn]{sn::rsn}.}

\item{moderators}{Atomic integer. The number of moderators to simulate for
each study. Make sure that the number of moderators to be simulated is at
least as large as the number of moderators referred to in the model
parameter. Internally, the matrix of moderators is referred to as \code{"x"}.
Defaults to 5.}

\item{distribution}{Atomic character. The distribution of the moderators.
Can be set to either \code{"normal"} or \code{"bernoulli"}. Defaults to \code{"normal"}.}

\item{model}{Expression. An expression to specify the model from which to
simulate the mean true effect size, mu. This formula may use the terms \code{"es"}
(referring to the es parameter of the call to simulate_smd), and \verb{"x\[, \]"}
(referring to the matrix of moderators, x). Thus, to specify that the mean
effect size, mu, is a function of the effect size and the first moderator,
one would pass the value \code{model = "es * x\[ , 1\]"}.
Defaults to \verb{"es * x\[ , 1\]"}.}
}
\value{
List of length 4. The "training" element of this list is a data.frame
with k_train rows. The columns are the variance of the effect size, vi; the
effect size, yi, and the moderators, X. The "testing" element of this list is
a data.frame with k_test rows. The columns are the effect size, yi, and the
moderators, X. The "housekeeping" element of this list is a data.frame with
k_train + k_test rows. The columns are n, the sample size n for each
simulated study; mu_i, the mean true effect size for each simulated study;
and theta_i, the true effect size for each simulated study.
}
\description{
This function simulates a meta-analytic dataset based on the random-effects
model. The simulated effect size is Hedges' G, an estimator of the
Standardized Mean Difference (Hedges, 1981; Li, Dusseldorp, & Meulman, 2017).
The functional form of the model can be specified, and moderators can be
either normally distributed or Bernoulli-distributed. See Van Lissa, in
preparation, for a detailed explanation of the simulation procedure.
}
\examples{
set.seed(8)
simulate_smd()
simulate_smd(k_train = 50, distribution = "bernoulli")
simulate_smd(distribution = "bernoulli", model = "es * x[ ,1] * x[ ,2]")
}
