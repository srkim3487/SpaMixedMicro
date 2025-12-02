library(spdep)
library(MASS)
library(Matrix)
library(phytools)
# library(VGAM)
library(INLA) 
library(phyloseq)
library(ANCOMBC)
# library(Rlab)
library(dplyr)
library(ade4)          # for dudi.pca
library(adespatial)    # for multispati.randtest
library(compositions)  # for CLR
library(fdrtool)

######################################################## feature selection
z_scores_fun <- function(marg) {
  mean_val <- inla.emarginal(function(x) x, marg)
  var_val  <- inla.emarginal(function(x) x^2, marg) - mean_val^2
  sd_val   <- sqrt(var_val)
  mean_val / sd_val
}

lfdr_fun <- function(pm_margs_all, delta = 0.2, sign_taxa){
  z_scores <- sapply(pm_margs_all, z_scores_fun)
  lfdr_res <- fdrtool(z_scores, statistic = "normal", plot = FALSE, cutoff.method="locfdr")
  lfdr_vals <- lfdr_res$lfdr
  
  sel_lfdr <- which(lfdr_vals < delta)
  
  TP <- length(intersect(sel_lfdr, sign_taxa))
  FP <- length(setdiff(sel_lfdr, sign_taxa))
  FN <- length(setdiff(sign_taxa, sel_lfdr))
  TN <- J - length(sign_taxa) - FP
  
  TPR <- TP / length(sign_taxa)
  FPR <- FP / (J - length(sign_taxa))
  PPV <- if((TP + FP) == 0) 0 else TP / (TP + FP)
  NPV <- if((TN + FN) == 0) 0 else TN / (TN + FN)
  
  return(list(sel_lfdr = sel_lfdr, 
              metric= c(TPR = TPR, FPR = FPR, PPV = PPV, NPV = NPV)))
}

average_metric_list <- function(fdr) {
  metrics <- c("TPR","FPR","PPV","NPV")
  M <- length(fdr)
  temp = t(sapply(1:M, function(m){(fdr[[m]]$metric)}))
  
  return(colMeans(temp))
}

ancombc_metric_fun <- function(res_iid, alpha, sign_taxa){
  temp = res_iid %>%
    filter(q_X_li < alpha)
  
  sel_lfdr <- sub("^T", "", temp$taxon)
  
  TP <- length(intersect(sel_lfdr, sign_taxa))
  FP <- length(setdiff(sel_lfdr, sign_taxa))
  FN <- length(setdiff(sign_taxa, sel_lfdr))
  TN <- J - length(sign_taxa) - FP
  
  TPR <- TP / length(sign_taxa)
  FPR <- FP / (J - length(sign_taxa))
  PPV <- if((TP + FP) == 0) 0 else TP / (TP + FP)
  NPV <- if((TN + FN) == 0) 0 else TN / (TN + FN)
  
  return(list(sel_lfdr = sel_lfdr, 
              metric= c(TPR = TPR, FPR = FPR, PPV = PPV, NPV = NPV)))
}



######################################################## fit
fit_fun <- function(data_gen, 
                    beta_prec1, beta_prec2,
                    u_prec1, u_prec2,
                    u_lambda1, u_lambda2,
                    v_prec1, v_prec2,
                    e_prec1, e_prec2,
                    e_lambda1, e_lambda2
                    ){
  
  X <- data_gen$X_li
  adj_mat <- data_gen$adj_mat
  micro_adj_mat <- data_gen$adj_mat_micro
  
  physeq <- data_gen$physeq
  micro_mat <- t(otu_table(physeq)) # n*J matrix
  micro_vec <- c(t(micro_mat))
  
  ######################################################## 
  N <- sum(indiv)
  Xmat <- kronecker(X*rep(1, N), diag(1, J))
  colnames(Xmat) <- paste("X", 1:ncol(Xmat), sep = "")
  
  data_df <- data.frame(micro_vec = micro_vec, 
                        int=factor(rep(1:J)),
                        beta_idx = factor(rep(1:J)),
                        Xmat,
                        Xdf = rep(X, each = J),
                        region=factor(rep(1:n, times = indiv * J)),
                        sampling_rnd=factor(rep(1:N, each = J)),
                        error = factor(rep(1:J, N)))
                        # offset = rep(offset_term, each = J))
  
  ######################################################## independent model
 res_ancom <- ancombc2(
   data = physeq,
   fix_formula = "X_li",            # covariate formula
   p_adj_method = "BH")

  ######################################################## spatial models
  regional_values <- as.factor(unique(levels(data_df$region)))
  micro_values <- as.factor(unique(levels(data_df$error)))
  
  formula <- paste("micro_vec ~ 0 + int +",
                        "+ f(beta_idx, Xdf, model='iid', hyper=list(prec = list(prior = 'pc.prec', param = c(", beta_prec1, ",", beta_prec2,"))))",
                        "+ f(region, values = regional_values, model = 'besagproper2', graph = adj_mat,
                                    hyper = list(prec = list(prior = 'pc.prec', param = c(", u_prec1, ",", u_prec2,")),
                                                 lambda = list(prior = 'logitbeta', param=c(", u_lambda1, ",", u_lambda2,"))))",
                   "+ f(sampling_rnd, model = 'iid', hyper=list(prec = list(prior = 'pc.prec', param = c(", v_prec1, ",", v_prec2,"))))",
                   "+ f(error, values = micro_values, model = 'besagproper2', graph=micro_adj_mat,
                                    hyper = list(prec = list(prior = 'pc.prec', param = c(", e_prec1, ",", e_prec2,")),
                                                 lambda = list(prior = 'logitbeta', param=c(", e_lambda1, ",", e_lambda2, "))))"
  )
  formula <- as.formula(formula)
  
  result <- suppressMessages(suppressWarnings(
    inla(formula,
         family = "zeroinflatedpoisson1", data = data_df, num.threads = 10,
         control.fixed = list(prec.intercept = 0.001),
         control.compute = list(config = TRUE, dic = TRUE, waic = TRUE))
  ))
 
  ###################################################################################### 
  pm_margs_all <- result$marginals.random$beta_idx
  
  return(list(data = data_gen, 
              sign_taxa = sign_taxa, betaM = betaM, 
              
              res_iid = res_ancom$res,
                        
              fixed_spa = result$summary.fixed, 
              betaA_spa = result$summary.random$beta_idx,
              rnd_spa = result$summary.hyperpar,
              
              pm_margs_all = pm_margs_all
  ))
}