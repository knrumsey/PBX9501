# Run RateStickAnalysis.R and RateStickAnalysis_SS_U_Ni.R first

# Get 50-fold cross validation measures for stainless steel and uranium
rmse <- function(zz) sqrt(mean(zz^2))
K <- 100 # folds
Nk <- 500/K
perm <- sample(1:500, 500, FALSE)
RMSE_U_S <- matrix(NA, 10, K)
for(k in 1:K){
  cat("Starting fold k = ", k, "\n")

  indx <- perm[((k-1)*Nk+1):(Nk*k)]
  ytrain_s <- y_m_id[[1]][-indx]
  ytest_s  <- y_m_id[[1]][indx]

  ytest_s  <- (ytest_s - mean(ytrain_s))/sd(ytrain_s)
  ytrain_s <- (ytrain_s - mean(ytrain_s))/sd(ytrain_s)

  ytrain_u <- y_m_id[[3]][-indx]
  ytest_u <- y_m_id[[3]][indx]

  ytest_u  <- (ytest_u - mean(ytrain_u))/sd(ytrain_u)
  ytrain_u <- (ytrain_u - mean(ytrain_u))/sd(ytrain_u)

  Xtrain <- X[-indx,]
  Xtest  <- X[indx,]

  Xtrain_u <- as.matrix(Xtrain)%*%matrix(eigen(C_u)$vectors[,1:2], nrow=6)
  Xtest_u  <- as.matrix(Xtest)%*%matrix(eigen(C_u)$vectors[,1:2], nrow=6)

  Xtrain_s <- as.matrix(Xtrain)%*%matrix(eigen(C_s)$vectors[,1:2], nrow=6)
  Xtest_s  <- as.matrix(Xtest)%*%matrix(eigen(C_s)$vectors[,1:2], nrow=6)

  Xtrain_su <- as.matrix(Xtrain)%*%matrix(eigen(V_su)$vectors[,c(1,6)], nrow=6)
  Xtest_su  <- as.matrix(Xtest)%*%matrix(eigen(V_su)$vectors[,c(1,6)], nrow=6)

  Xtrain_zahm <- as.matrix(Xtrain)%*%matrix(eigen(C_u+C_s)$vectors[,1:2], nrow=6)
  Xtest_zahm  <- as.matrix(Xtest)%*%matrix(eigen(C_u+C_s)$vectors[,1:2], nrow=6)


  # Fit models
  mod_s <- bass(Xtrain, ytrain_s, verbose=FALSE)
  mod_u <- bass(Xtrain, ytrain_u, verbose=FALSE)

  mod_us <- bass(Xtrain_u, ytrain_s, verbose=FALSE)
  mod_uu <- bass(Xtrain_u, ytrain_u, verbose=FALSE)
  cat("\t --")

  mod_ss <- bass(Xtrain_s, ytrain_s, verbose=FALSE)
  mod_su <- bass(Xtrain_s, ytrain_u, verbose=FALSE)
  cat("--")

  mod_sus <- bass(Xtrain_su, ytrain_s, verbose=FALSE)
  mod_suu <- bass(Xtrain_su, ytrain_u, verbose=FALSE)
  cat("--\n")

  mod_zs <- bass(Xtrain_zahm, ytrain_s, verbose=FALSE)
  mod_zu <- bass(Xtrain_zahm, ytrain_u, verbose=FALSE)

  # Get predictions
  RMSE_U_S[1,k] <- rmse(apply(predict(mod_us, Xtest_u),     2, mean) - ytest_s)
  RMSE_U_S[2,k] <- rmse(apply(predict(mod_uu, Xtest_u),     2, mean) - ytest_u)
  RMSE_U_S[3,k] <- rmse(apply(predict(mod_ss, Xtest_s),     2, mean) - ytest_s)
  RMSE_U_S[4,k] <- rmse(apply(predict(mod_su, Xtest_s),     2, mean) - ytest_u)
  RMSE_U_S[5,k] <- rmse(apply(predict(mod_sus, Xtest_su),   2, mean) - ytest_s)
  RMSE_U_S[6,k] <- rmse(apply(predict(mod_suu, Xtest_su),   2, mean) - ytest_u)
  RMSE_U_S[7,k] <- rmse(apply(predict(mod_zs, Xtest_zahm),  2, mean) - ytest_s)
  RMSE_U_S[8,k] <- rmse(apply(predict(mod_zu, Xtest_zahm),  2, mean) - ytest_u)

  #Full data
  RMSE_U_S[9,k] <- rmse(apply(predict(mod_s, Xtest), 2, mean) - ytest_s)
  RMSE_U_S[10,k] <- rmse(apply(predict(mod_u, Xtest), 2, mean) - ytest_u)
}
tmp <- apply(RMSE_U_S,1,mean)
TAB <- matrix(NA, nrow=2, ncol=5)
TAB[,1] <- tmp[9:10]
TAB[,2] <- tmp[3:4]
TAB[,3] <- tmp[1:2]
TAB[,4] <- tmp[5:6]
TAB[,5] <- tmp[7:8]
rownames(TAB) <- c("Stainless Steel 304", "Uranium")
colnames(TAB) <- c("Full Data", "Css", "Cu", "Vus", "Cu+Cs")
stargazer::stargazer(TAB)


