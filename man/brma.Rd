% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/brma.R
\name{brma}
\alias{brma}
\alias{brma.formula}
\alias{brma.default}
\title{Conduct Bayesian Regularized Meta-Analysis}
\usage{
brma(x, ...)

\method{brma}{formula}(
  formula,
  data,
  vi = "vi",
  study = NULL,
  method = "hs",
  standardize = TRUE,
  prior = switch(method, lasso = c(df = 1, scale = 1), hs = c(df = 1, df_global = 1,
    df_slab = 4, scale_global = 1, scale_slab = 2, relevant_pars = NULL)),
  mute_stan = TRUE,
  ...
)

\method{brma}{default}(
  x,
  y,
  vi,
  study = NULL,
  method = "hs",
  standardize,
  prior,
  mute_stan = TRUE,
  intercept,
  ...
)
}
\arguments{
\item{x}{An k x m numeric matrix, where k is the number of effect sizes and m
is the number of moderators.}

\item{...}{Additional arguments passed on to \code{\link[rstan:stanmodel-method-sampling]{rstan::sampling()}}.
Use this, e.g., to override default arguments of that function.}

\item{formula}{An object of class \code{formula} (or one that can be coerced to
that class), see \code{\link[stats]{lm}}.}

\item{data}{Either a \code{data.frame} containing the variables in the model,
see \code{\link[stats]{lm}}, or a \code{list} of multiple imputed \code{data.frame}s,
or an object returned by \code{\link[mice]{mice}}.}

\item{vi}{Character. Name of the column in the \code{data} that
contains the variances of the effect sizes. This column will be removed from
the data prior to analysis. Defaults to \code{"vi"}.}

\item{study}{Character. Name of the column in the
\code{data} that contains the study id. Use this when the data includes
multiple effect sizes per study. This column can be a vector of integers, or
a factor. This column will be removed from the data prior to analysis.
See \code{Details} for more information about analyzing dependent data.}

\item{method}{Character, indicating the type of regularizing prior to use.
Supports one of \code{c("hs", "lasso")}, see Details. Defaults to
\code{"hs"}.}

\item{standardize}{Either a logical argument or a list. If \code{standardize} is
logical, it controls whether all predictors are standardized prior to
analysis or not. Parameter estimates are restored to the predictors' original
scale. Alternatively, users can provide a list to \code{standardize} to gain
more control over the standardization process. In this case, it is assumed
that the standardization has already taken place. This list must have two
elements: \code{list(center = c(mean(X1) , mean(X2), mean(X...)), scale = c(sd(X1), sd(X2), sd(X...)))}. It is used
only to restore parameter estimates to the original scale of the predictors.
This is useful, e.g., to standardize continuous and dichotomous variables
separately. Defaults to \code{TRUE}, which is recommended so that shrinking
affects all parameters similarly.}

\item{prior}{Numeric vector, specifying the prior to use. Note that the
different \code{method}s require this vector to contain specific named
elements.}

\item{mute_stan}{Logical, indicating whether mute all 'Stan' output or not.}

\item{y}{A numeric vector of k effect sizes.}

\item{intercept}{Logical, indicating whether or not an intercept should be included
in the model.}
}
\value{
A \code{list} object of class \code{brma}, with the following structure:

\if{html}{\out{<div class="sourceCode">}}\preformatted{list(
  fit          # An object of class stanfit, for compatibility with rstan
  coefficients # A numeric matrix with parameter estimates; these are
               # interpreted as regression coefficients, except tau2 and tau,
               # which are interpreted as the residual variance and standard
               # deviation, respectively.
  formula      # The formula used to estimate the model
  terms        # The predictor terms in the formula
  X            # Numeric matrix of moderator variables
  Y            # Numeric vector with effect sizes
  vi           # Numeric vector with effect size variances
  tau2         # Numeric, estimated tau2
  R2           # Numeric, estimated heterogeneity explained by the moderators
  k            # Numeric, number of effect sizes
  study        # Numeric vector with study id numbers
)
}\if{html}{\out{</div>}}
}
\description{
This function conducts Bayesian regularized meta-regression (Van Lissa & Van
Erp, 2021). It uses the \code{stan} function
\link[rstan:stanmodel-method-sampling]{rstan::sampling} to fit the model. A lasso or horseshoe prior is used to
shrink the regression coefficients of irrelevant moderators towards zero.
See Details.
}
\details{
The Bayesian regularized meta-analysis algorithm (Van Lissa & Van
Erp, 2021) penalizes meta-regression coefficients either via the
lasso prior (Park & Casella, 2008) or the regularized horseshoe prior
(Piironen & Vehtari, 2017).
\describe{
\item{lasso}{ The Bayesian equivalent of the lasso penalty is obtained when
placing independent Laplace (i.e., double exponential) priors on the
regression coefficients centered around zero. The scale of the Laplace
priors is determined by a global scale parameter \code{scale}, which
defaults to 1 and an inverse-tuning parameter \eqn{\frac{1}{\lambda}}
which is given a chi-square prior governed by a degrees of freedom
parameter \code{df} (defaults to 1). If \code{standardize = TRUE},
shrinkage will
affect all coefficients equally and it is not necessary to adapt the
\code{scale} parameter. Increasing the \code{df} parameter will allow
larger values for the inverse-tuning parameter, leading to less shrinkage.}
\item{hs}{ One issue with the lasso prior is that it has relatively light
tails. As a result, not only does the lasso have the desirable behavior of
pulling small coefficients to zero, it also results in too much shrinkage
of large coefficients. An alternative prior that improves upon this
shrinkage pattern is the horseshoe prior (Carvalho, Polson & Scott, 2010).
The horseshoe prior has an infinitely large spike at zero, thereby pulling
small coefficients toward zero but in addition has fat tails, which allow
substantial coefficients to escape the shrinkage. The regularized horseshoe
is an extension of the horseshoe prior that allows the inclusion of prior
information regarding the number of relevant predictors and can
be more numerically stable in certain cases (Piironen & Vehtari, 2017).
The regularized horseshoe has a global shrinkage parameter that influences
all coefficients similarly and local shrinkage parameters that enable
flexible shrinkage patterns for each coefficient separately. The local
shrinkage parameters are given a Student's t prior with a default \code{df}
parameter of 1. Larger values for \code{df} result in lighter tails and
a prior that is no longer strictly a horseshoe prior. However, increasing
\code{df} slightly might be necessary to avoid divergent transitions in
Stan (see also \url{https://mc-stan.org/misc/warnings.html}). Similarly,
the degrees of freedom for the Student's t prior on the global shrinkage
parameter \code{df_global} can be increased from the default of 1 to, for
example, 3 if divergent transitions occur although the resulting
prior is then strictly no longer a horseshoe. The scale for the Student's t
prior on the global shrinkage parameter \code{scale_global} defaults to 1
and can be decreased to achieve more shrinkage. Moreover, if prior
information regarding the number of relevant moderators is available, it is
recommended to include this information via the \code{relevant_pars}
argument by setting it to the expected number of relevant moderators. When
\code{relevant_pars} is specified, \code{scale_global} is ignored and
instead based on the available prior information. Contrary to the horseshoe
prior, the regularized horseshoe applies additional regularization on large
coefficients which is governed by a Student's t prior with a
\code{scale_slab} defaulting to 2 and \code{df_slab} defaulting to 4.
This additional regularization ensures at least some shrinkage of large
coefficients to avoid any sampling problems.}
}
}
\examples{
data("curry")
df <- curry[c(1:5, 50:55), c("d", "vi", "sex", "age", "donorcode")]
suppressWarnings({res <- brma(d~., data = df, iter = 10)})
}
\references{
Van Lissa, C. J., van Erp, S., & Clapper, E. B. (2023). Selecting relevant
moderators with Bayesian regularized meta-regression. Research Synthesis
Methods. \doi{10.31234/osf.io/6phs5}

Park, T., & Casella, G. (2008). The Bayesian Lasso. Journal of the American
Statistical Association, 103(482), 681–686. \doi{10.1198/016214508000000337}

Carvalho, C. M., Polson, N. G., & Scott, J. G. (2010). The horseshoe
estimator for sparse signals. Biometrika, 97(2), 465–480.
\doi{10.1093/biomet/asq017}

Piironen, J., & Vehtari, A. (2017). Sparsity information and regularization
in the horseshoe and other shrinkage priors. Electronic Journal of
Statistics, 11(2). \url{https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-11/issue-2/Sparsity-information-and-regularization-in-the-horseshoe-and-other-shrinkage/10.1214/17-EJS1337SI.pdf}
}
