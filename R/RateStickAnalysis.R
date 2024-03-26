library(concordance) # devtools::install_github("knrumsey/concordance")
library(BASS)
library(stargazer)
library(fields)
tr <- function(M) sum(diag(M))

jacket_materials <- jacket_names <- c("copper", "tungsten", "ss304", "gold", "gold_5cu",
                      "uranium_5mo", "nickel", "uranium_075ti", "al6061",
                      "uranium", "ss250", "tin", "ss4340", "al7075")
fnames <- unlist(lapply(1:14, function(zz) paste0(jacket_materials[zz],".csv")))
p <- 6 # Number of Variables

########################################################
#             FIT BASS MODELS AND ESTIMATE C
########################################################
mod <- C <-  list()
mc_ind <- seq(100, 1000, length=10)
cnt <- 1
for(j in 1:length(fnames)){
  tmp <- read.csv(paste0("data/", fnames[j]))
  X <- tmp[, 12:17]
  ord <- order(X[,1])
  X <- X[ord,]
  for(i in 1:ncol(X)) X[,i] <- BASS:::scale.range(X[,i])
  y <- tmp[ord,11]
  mod[[j]] <- bass(X, y, h1=1, h2=0.1)
  plot(mod[[j]])
  for(k in seq_along(mc_ind)){
    C[[cnt]] <- C_bass(mod[[j]], mcmc.use=mc_ind[k])
    cnt <- cnt + 1
  }
}

Cij <- list()
hold_i <- hold_j <- rep(NA, choose(140, 2))
CONC <- matrix(NA, nrow=length(mod)*length(mc_ind), ncol=length(mod)*length(mc_ind))
cnt <- 1
for(i in 2:140){
  j1 <- ceiling(i/10)
  k1 <- (i %% 10) + 1
  for(j in 1:(i-1)){
    j2 <- ceiling(j/10)
    k2 <- (j %% 10) + 1
    Cij[[cnt]] <- Cfg_bass(mod[[j1]], mod[[j2]], mcmc.use=cbind(k1, k2))
    CONC[i,j] <- CONC[j,i] <- tr(Cij[[cnt]])/sqrt(tr(C[[i]])*tr(C[[j]]))
    hold_i[cnt] <- j1
    hold_j[cnt] <- j2
    cnt <- cnt + 1
    if(cnt %% 250 == 0) cat(cnt, "\n")
  }
}
CONC <- pmin(CONC, 0.99999) # numerical stability
diag(CONC) <- 1

CONC_mat <- array(NA, dim=c(14, 14, 100))
CONC_mat_mean <- array(NA, dim=c(14, 14))
CONC_mat_sd <- array(NA, dim=c(14, 14))
for(i in 1:14){
  for(j in 1:14){
    ind_i <- ((i-1)*10+1):(10*i)
    ind_j <- ((j-1)*10+1):(10*j)
    CONC_mat[i,j,] <- conc_post <- as.numeric(CONC[ind_i, ind_j])
    CONC_mat_mean[i,j] <- mean(conc_post)
    CONC_mat_sd[i,j] <- sd(conc_post)
  }
}

########################################################
#         MAKE HEAT MAPS FOR CONCORDANCE VALUES
########################################################
png("figs/image_conc_mean.png", units="in", height=5, width=5.5, res=300)
par(mar=c(5.1, 4.1, 4.1, 2.1) + c(1.5, 2.2, 0, 0))
image.plot(CONC_mat_mean, xaxt='n', yaxt='n', col=hcl.colors(50, "YlOrRd", rev=TRUE),
           main="Concordance (Posterior mean)")
axis(1, seq(0, 1, length=14), jacket_names, las=2)
axis(2, seq(0, 1, length=14), jacket_names, las=1)
dev.off()

png("figs/image_conc_sd.png", units="in", height=5, width=5.5, res=300)
par(mar=c(5.1, 4.1, 4.1, 2.1) + c(1.5, 2.2, 0, 0))
image.plot(CONC_mat_sd, xaxt='n', yaxt='n', col=hcl.colors(50, "YlOrRd", rev=TRUE),
           main="Concordance (Posterior sd)")
axis(1, seq(0, 1, length=14), jacket_names, las=2)
axis(2, seq(0, 1, length=14), jacket_names, las=1)
dev.off()

png("figs/image_conc_samples.png", units="in", height=5, width=5.5, res=300)
par(mar=c(5.1, 4.1, 4.1, 2.1) + c(1.5, 2.2, 0, 0))
image.plot(CONC, xaxt='n', yaxt='n', col=hcl.colors(50, "YlOrRd", rev=TRUE),
           main="Concordance (Posterior samples)")
axis(1, seq(0.025, 0.975, length=14), jacket_names, las=2)
axis(2, seq(0.025, 0.975, length=14), jacket_names, las=1)
dev.off()

########################################################
#             MAKE VORONOI DIAGRAM
########################################################
library(MASS)
library(tripack)
library(RColorBrewer)
bob <- c(brewer.pal(12, "Paired"), "grey70", "firebrick4")
bob[11] <- "gold"
bob <- bob[c(7,2,12,11,8,5,3,9,1,10,6,14,4,13)] #Customize colors for better presentation
DISC <- sqrt(1 - CONC)

x <- MASS::isoMDS(DISC)$points
z <- matrix(NA, nrow=14, ncol=10)
for(i in 1:14){
  indx <- ((i-1)*10+1):(10*i)
  z[i,] <- apply(x[indx,], 2, mean)
}
vm <- tripack::voronoi.mosaic(z[,1], z[,2])

png(filename = "figs/voronoi_post.png", height=5, width=8, res=400, units="in")
par(mar=c(.1, .1, .1, .1))
plot(vm, main="", sub="", do.points=FALSE, lwd=2, col='gray')
points(x, pch=16, col=adjustcolor(bob[rep(1:14, each=10)], alpha.f=0.4))
points(z, pch=21, bg=bob, cex=1.4)
n <- length(vm$x)
cnt <- 0
for(i in 1:n){
  if (vm$node[i]) {
    tns <- sort(c(vm$n1[i], vm$n2[i], vm$n3[i]))
    for(j in 1:3){
      if(tns[j] < 0)
        lines(c(vm$x[i], vm$dummy.x[-tns[j]]), c(vm$y[i], vm$dummy.y[-tns[j]]), lty = 1, lwd=2, col='gray')
    }
  }
}
leave_blank <- c(3, 9, 11, 13, 14)
eps <- 0
# ALUMINUM ARROWS
nnn <- 2
arrows(z[leave_blank[2],1], z[leave_blank[2],2], 3.75, 0,
       col=bob[leave_blank[2]], angle=18, length=0.15, lwd=1, code=1, lty=2)
text(3.75, 0 + eps, labels=jacket_names[leave_blank[2]], cex=0.9, font=4)
#text(3.75, 0 + eps, labels=jacket_names[leave_blank[2]], cex=0.78, font=4, col=bob[leave_blank[nnn]])

nnn <- 5
arrows(z[leave_blank[5],1], z[leave_blank[5],2], 3.75, -0.35,
       col=bob[leave_blank[5]], angle=18, length=0.15, lwd=1, code=1, lty=2)
text(3.75, -0.35 + eps, labels=jacket_names[leave_blank[5]], cex=0.9, font=4)
#text(3.75, -0.55 + eps, labels=jacket_names[leave_blank[5]], cex=0.78, font=4, col=bob[leave_blank[nnn]])

# SS ARROWS
nnn <- 4
arrows(z[leave_blank[nnn],1], z[leave_blank[nnn],2], 3.75, -.65,
       col=bob[leave_blank[nnn]], angle=18, length=0.15, lwd=1, code=1, lty=2)
text(3.75, -.65 + eps, labels=jacket_names[leave_blank[nnn]], cex=0.9, font=4)
#text(3.75, -1.3 + eps, labels=jacket_names[leave_blank[nnn]], cex=0.78, font=4, col=bob[leave_blank[nnn]])

nnn <- 1
arrows(z[leave_blank[nnn],1], z[leave_blank[nnn],2], -2.4, 1.8,
       col=bob[leave_blank[nnn]], angle=18, length=0.15, lwd=1, code=1, lty=2)
text(-2.4, 1.8 - eps, labels=jacket_names[leave_blank[nnn]], cex=0.9, font=4)
#text(-3, 1.8 - eps, labels=jacket_names[leave_blank[nnn]], cex=0.78, font=4, col=bob[leave_blank[nnn]])

nnn <- 3
arrows(z[leave_blank[nnn],1], z[leave_blank[nnn],2], -1.6, 1.8,
       col=bob[leave_blank[nnn]], angle=18, length=0.15, lwd=1, code=1, lty=2)
text(-1.5, 1.8 - eps, labels=jacket_names[leave_blank[nnn]], cex=0.9, font=4)
#text(-1.5, 1.8 - eps, labels=jacket_names[leave_blank[nnn]], cex=0.78, font=4, col=bob[leave_blank[nnn]])
# Add remaining text
eps <- matrix(-.15, nrow=9, ncol=2)
eps[1,] <- c(0.03, -.15)
text(z[-leave_blank,1:2] + eps, labels=jacket_names[-leave_blank], cex=0.9, font=4)
#text(z[-leave_blank,] + eps, labels=jacket_names[-leave_blank], cex=0.78, font=4, col=bob[-leave_blank])
dev.off()
par(mar=c(5.1, 4.1, 4.1, 2.1))


########################################################
#        MAKE TABLES OF CONCORDANCE VALUES
########################################################
library(stargazer)
rownames(CONC_mat_mean) <- jacket_names
colnames(CONC_mat_mean) <- jacket_names
rownames(CONC_mat_sd) <- jacket_names
colnames(CONC_mat_sd) <- jacket_names

stargazer(CONC_mat_mean[,1:7], summary=FALSE)
stargazer(CONC_mat_mean[,8:14], summary=FALSE)
stargazer(CONC_mat_sd[,1:7], summary=FALSE)
stargazer(CONC_mat_sd[,8:14], summary=FALSE)

#Asses model fits
K <- 100 # folds
Nk <- 500/K
RMSE <- matrix(NA, nrow=14, ncol=K)
for(i in 1:14){
  tmp <- read.csv(paste0("data/", fnames[i]))
  X <- tmp[, 12:17]
  ord <- order(X[,1])
  X <- X[ord,]
  for(ii in 1:ncol(X)) X[,ii] <- BASS:::scale.range(X[,ii])
  y <- tmp[ord,11]
  y <- (y-mean(y))/sd(y)

  #Cross validation

  perm <- sample(1:500, 500, FALSE)
  for(k in 1:K){
    indx <- perm[((k-1)*Nk+1):(Nk*k)]
    ytrain <- y[-indx]
    ytest  <- y[indx]
    Xtrain <- X[-indx,]
    Xtest  <- X[indx,]
    mod_curr <- bass(Xtrain, ytrain, h1=1, h2=0.1)
    yhat <- apply(predict(mod_curr, Xtest), 2, mean)
    RMSE[i,k] <- sqrt(mean((yhat-ytest)^2))
  }
  cat("\n\n\n\n\n\n", i, "\n\n\n\n\n\n")
}
cv_rmse = (round(apply(RMSE, 1, mean), 4))
range(cv_rmse)
jacket_names[order(cv_rmse)]

