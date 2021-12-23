plot_prior <- function(method = c("hs", "lasso"),
                       prior = switch(method,
                                      "lasso" = c(df = 1, scale = 1),
                                      "hs" = c(df = 1, df_global = 1, df_slab = 4, scale_global = 1, scale_slab = 1, par_ratio = NULL)),
                       iter = 1000){
  smpls <- suppressWarnings(sampling(object = pema:::stanmodels[[c("lasso_prior", "hs_prior")[(method[1] == "hs")+1]]], data = as.list(prior), chains = 1, iter = 10000, warmup = 0, show_messages = FALSE, verbose = FALSE, refresh = 0))
  plot(density(smpls@sim$samples[[1]]$b), main = c("Lasso prior", "Horseshoe prior")[(method == "hs")+1],
       xlab = paste0("Samples: ", iter, ", ", paste0(names(prior), " = ", prior, collapse = ", ")),
       xlim = c(-5, 5))
}
