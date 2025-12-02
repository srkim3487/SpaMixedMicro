######################################################## 
source("utils.R")
load("sim_data_example.RData")
######################################################## 
adj_mat_micro <- data_gen$adj_mat_micro
sign_taxa <- data_gen$sign_taxa
beta0 <- data_gen$beta0
betaM <- data_gen$betaM
indiv <- data_gen$indiv
J <- length(betaM)
n <- dim(data_gen$data_gen$adj_mat)[1]

data <- data_gen$data_gen
######################################################## Fit
res <- fit_fun(data, 
         beta_prec1=0.5, beta_prec2=0.01,
         u_prec1 = 2, u_prec2 = 0.05,
         u_lambda1 = 1, u_lambda2 = 1,
         v_prec1 = 2, v_prec2 = 0.05,
         e_prec1 = 2, e_prec2 = 0.05,
         e_lambda1 = 2, e_lambda2 = 2)

######################################################## feature selection
## ANCOM-BC2
ancombc_metric_fun(res$res_iid, alpha = 0.1, sign_taxa)

## spatial model
lfdr_fun(res$pm_margs_all, delta = 0.2, sign_taxa)

