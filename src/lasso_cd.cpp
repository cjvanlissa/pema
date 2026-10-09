#include <Rcpp.h>
#include <cmath>
#include <algorithm>

using namespace Rcpp;


// Soft-thresholding operator
inline double soft_threshold(double z, double lambda) {
  if (z > lambda) {
    return z - lambda;
  } else if (z < -lambda) {
    return z + lambda;
  }
  return 0.0;
}


// [[Rcpp::export]]
List lasso_cd_cpp(
    const NumericMatrix& X,
    const NumericVector& y,
    const NumericVector& w,
    double lambda,
    bool intercept = true,
    Nullable<NumericVector> beta_start = R_NilValue,
    double tol = 1e-8,
    int maxit = 10000
) {

  const int n = X.nrow();
  const int p = X.ncol();

  // Input validation
  if (y.size() != n || w.size() != n) {
    stop("Incompatible dimensions of X, y, and w.");
  }

  if (p < 1) {
    stop("X must contain at least one column.");
  }

  if (!std::isfinite(lambda) || lambda < 0.0) {
    stop("lambda must be finite and nonnegative.");
  }

  if (!std::isfinite(tol) || tol <= 0.0) {
    stop("tol must be finite and positive.");
  }

  if (maxit < 1) {
    stop("maxit must be at least 1.");
  }

  // Initialize coefficients
  NumericVector beta(p);

  if (beta_start.isNotNull()) {
    NumericVector initial(beta_start);

    if (initial.size() != p) {
      stop("beta_start must have length ncol(X).");
    }

    std::copy(initial.begin(), initial.end(), beta.begin());
  }

  // Validate data and weights
  for (int i = 0; i < n; ++i) {
    if (!std::isfinite(y[i]) ||
        !std::isfinite(w[i]) ||
        w[i] <= 0.0) {
      stop("y must be finite and w must be positive and finite.");
    }
  }

  for (int j = 0; j < p; ++j) {
    if (!std::isfinite(beta[j])) {
      stop("beta_start must contain only finite values.");
    }

    for (int i = 0; i < n; ++i) {
      if (!std::isfinite(X(i, j))) {
        stop("X must contain only finite values.");
      }
    }
  }

  // Initialize residuals: r = y - X * beta
  NumericVector r = clone(y);

  for (int j = 0; j < p; ++j) {
    if (beta[j] != 0.0) {
      for (int i = 0; i < n; ++i) {
        r[i] -= X(i, j) * beta[j];
      }
    }
  }

  // Weighted squared column norms:
  // z[j] = sum_i w[i] * X(i,j)^2
  NumericVector z(p);

  for (int j = 0; j < p; ++j) {
    double zj = 0.0;

    for (int i = 0; i < n; ++i) {
      const double xij = X(i, j);
      zj += w[i] * xij * xij;
    }

    z[j] = zj;
  }

  if (intercept && z[0] <= 0.0) {
    stop("The intercept column must have positive weighted norm.");
  }

  bool converged = false;
  int iterations = 0;

  // Coordinate descent
  for (int iter = 0; iter < maxit; ++iter) {

    iterations = iter + 1;

    NumericVector beta_old = clone(beta);

    // --------------------------------------------------
    // Update unpenalized intercept (first column)
    // --------------------------------------------------
    if (intercept) {

      const int j = 0;

      // Add back current coefficient contribution
      for (int i = 0; i < n; ++i) {
        r[i] += X(i, j) * beta[j];
      }

      double rho = 0.0;

      for (int i = 0; i < n; ++i) {
        rho += w[i] * X(i, j) * r[i];
      }

      beta[j] = rho / z[j];

      // Remove updated contribution
      for (int i = 0; i < n; ++i) {
        r[i] -= X(i, j) * beta[j];
      }
    }

    // --------------------------------------------------
    // Update penalized moderator coefficients
    // --------------------------------------------------
    const int first_penalized = intercept ? 1 : 0;

    for (int j = first_penalized; j < p; ++j) {

      // Add back current coefficient contribution
      for (int i = 0; i < n; ++i) {
        r[i] += X(i, j) * beta[j];
      }

      // Weighted partial correlation
      double rho = 0.0;

      for (int i = 0; i < n; ++i) {
        rho += w[i] * X(i, j) * r[i];
      }

      // Soft-thresholding update
      if (z[j] > 0.0) {
        beta[j] = soft_threshold(rho, lambda) / z[j];
      } else {
        beta[j] = 0.0;
      }

      // Update residuals
      for (int i = 0; i < n; ++i) {
        r[i] -= X(i, j) * beta[j];
      }
    }

    // --------------------------------------------------
    // Convergence criterion
    // --------------------------------------------------
    double max_change = 0.0;
    double max_old = 0.0;

    for (int j = 0; j < p; ++j) {
      max_change = std::max(
        max_change,
        std::abs(beta[j] - beta_old[j])
      );

      max_old = std::max(
        max_old,
        std::abs(beta_old[j])
      );
    }

    if (max_change < tol * (1.0 + max_old)) {
      converged = true;
      break;
    }

    // Allow user interruption
    if (iter % 100 == 0) {
      Rcpp::checkUserInterrupt();
    }
  }

  // Preserve coefficient names if supplied
  if (X.hasAttribute("dimnames")) {
    List dimnames = X.attr("dimnames");

    if (!Rf_isNull(dimnames[1])) {
      beta.attr("names") = dimnames[1];
    }
  }

  return List::create(
    Named("beta") = beta,
    Named("residuals") = r,
    Named("iterations") = iterations,
    Named("converged") = converged
  );
}
