# NOTE: Run RateStickAnalysis.R first

########################################################
#        GET CO-SENSITIVITIES
########################################################

# Co-Sensitivity Analysis
m_id <- c(3, 7, 10) # ss304, nickel and uranium
m_nms <- jacket_names[m_id]

# Get Constantine matrices
j <- m_id[1]
C_s <- matrix(0, p, p)
for(jj in ((j-1)*10+1):(10*j)){
  C_s <- C_s + C[[jj]]/10
}

j <- m_id[2]
C_n <- matrix(0, p, p)
for(jj in ((j-1)*10+1):(10*j)){
  C_n <- C_n + C[[jj]]/10
}

j <- m_id[3]
C_u <- matrix(0, p, p)
for(jj in ((j-1)*10+1):(10*j)){
  C_u <- C_u + C[[jj]]/10
}

# Get Symmetrized Matrices
i <- max(m_id[1:2])
j <- min(m_id[1:2])
indx <- intersect(which(hold_i == i), which(hold_j == j))
C_sn <- matrix(0, p, p)
for(k in indx){
  C_sn <- C_sn + Cij[[k]]/100
}
V_sn <- (C_sn + t(C_sn))/2

i <- max(m_id[2:3])
j <- min(m_id[2:3])
indx <- intersect(which(hold_i == i), which(hold_j == j))
C_nu <- matrix(0, p, p)
for(k in indx){
  C_nu <- C_nu + Cij[[k]]/100
}
V_nu <- (C_nu + t(C_nu))/2

i <- max(m_id[c(1,3)])
j <- min(m_id[c(1,3)])
indx <- intersect(which(hold_i == i), which(hold_j == j))
C_su <- matrix(0, p, p)
for(k in indx){
  C_su <- C_su + Cij[[k]]/100
}
V_su <- (C_su + t(C_su))/2


########################################################
#                    MAKE TABLE 3
########################################################
pi_sn <- eigen(V_sn)$values/sqrt(tr(C_s)*tr(C_n))
sum(pi_sn) # concordance

pi_nu <- eigen(V_nu)$values/sqrt(tr(C_u)*tr(C_n))
sum(pi_nu)

pi_su <- eigen(V_su)$values/sqrt(tr(C_u)*tr(C_s))
sum(pi_su)


########################################################
#                    MAKE TABLE 4
########################################################
A <- matrix(NA, nrow=p, ncol=6)
A[,1] <- act_scores(C_s)
A[,2] <- act_scores(C_n)
A[,3] <- act_scores(C_u)
A[,4] <- coact_scores(V_sn)
A[,5] <- coact_scores(V_nu)
A[,6] <- coact_scores(V_su)
A <- sign(A)*sqrt(abs(A))
rownames(A) <- c("r0", "a", "b", "r1", "r2", "w")
colnames(A) <- c("Cs", "Cn", "Cu", "V_sn", "V_nu", "V_su")
stargazer(t(1000*A))


########################################################
#                    MAKE FIGURE 7
########################################################
png("figs/activity1.png", units="in", height=5, width=5, res=300)
plot(A[,3]/A[,1], type='o', pch=21, bg="orange",
     lwd=2, lty=1, col="orange", cex=1.3,
     xaxt="n", xlab="", ylab="Relative Activity")
#points(A[,3]/A[,1], cex=1.3, lwd=2)
lines(A[,2]/A[,1], lwd=2, col="dodgerblue")
points(A[,2]/A[,1], pch=15, col="dodgerblue", cex=1.3, lwd=2)
abline(h=1, lty=3)
axis(1, 1:6, c(expression(paste(rho[0])), "A", "B", expression(paste("r"[1])), expression(paste("r"[2])), expression(omega)))
legend("bottomleft",
       c(expression(paste(over(alpha["U"],alpha["SS"]))),
         expression(paste(over(alpha["Ni"],alpha["SS"])))),
       lwd=2, col=c("orange", "dodgerblue"), pch=15:16, cex=1.2, horiz=TRUE, bty='y', inset=c(0,0))
dev.off()

png("figs/activity2.png", units="in", height=5, width=5, res=300)
xx1 <- A[,1]/A[,4]
xx2 <- A[,2]/A[,4]
plot(xx1, type='o', pch=21, bg="orange",
     lwd=2, lty=1, col="orange", cex=1.3,
     xaxt="n", xlab="", ylab="Relative Activity", ylim=c(0.65, 1.3))#range(c(xx1, xx2)))
#points(A[,3]/A[,1], cex=1.3, lwd=2)
lines(xx2, lwd=2, col="dodgerblue")
points(xx2, pch=15, col="dodgerblue", cex=1.3, lwd=2)
abline(h=1, lty=3)
axis(1, 1:6, c(expression(paste(rho[0])), "A", "B", expression(paste("r"[1])), expression(paste("r"[2])), expression(omega)))
legend("bottomleft",
       c(expression(paste(over(alpha["SS"],alpha["SS,Ni"]))),
         expression(paste(over(alpha["Ni"],alpha["SS,Ni"])))),
       lwd=2, col=c("orange", "dodgerblue"), pch=15:16, cex=1.2, horiz=TRUE, bty='y', inset=c(0,0))
dev.off()

png("figs/activity3.png", units="in", height=5, width=5, res=300)
xx1 <- A[,2]/A[,5]
xx2 <- A[,3]/A[,5]
plot(xx1, type='o', pch=21, bg="orange",
     lwd=2, lty=1, col="orange", cex=1.3,
     xaxt="n", xlab="", ylab="Relative Activity", ylim=range(c(xx1, xx2)))
#points(A[,3]/A[,1], cex=1.3, lwd=2)
lines(xx2, lwd=2, col="dodgerblue")
points(xx2, pch=15, col="dodgerblue", cex=1.3, lwd=2)
abline(h=1, lty=3)
axis(1, 1:6, c(expression(paste(rho[0])), "A", "B", expression(paste("r"[1])), expression(paste("r"[2])), expression(omega)))
legend("bottomleft",
       c(expression(paste(over(alpha["Ni"],alpha["Ni,U"]))),
         expression(paste(over(alpha["U"],alpha["Ni,U"])))),
       lwd=2, col=c("orange", "dodgerblue"), pch=15:16, cex=1.2, horiz=TRUE, bty='y', inset=c(0,0))
dev.off()


png("figs/activity4.png", units="in", height=5, width=5, res=300)
xx1 <- A[,1]/A[,6]
xx2 <- A[,3]/A[,6]
plot(xx1, type='o', pch=21, bg="orange",
     lwd=2, lty=1, col="orange", cex=1.3,
     xaxt="n", xlab="", ylab="Relative Activity", ylim=c(0.65, 1.3))#range(c(xx1, xx2)))
#points(A[,3]/A[,1], cex=1.3, lwd=2)
lines(xx2, lwd=2, col="dodgerblue")
points(xx2, pch=15, col="dodgerblue", cex=1.3, lwd=2)
abline(h=1, lty=3)
axis(1, 1:6, c(expression(paste(rho[0])), "A", "B", expression(paste("r"[1])), expression(paste("r"[2])), expression(omega)))
legend("bottomleft",
       c(expression(paste(over(alpha["SS"],alpha["SS,U"]))),
         expression(paste(over(alpha["U"],alpha["SS,U"])))),
       lwd=2, col=c("orange", "dodgerblue"), pch=15:16, cex=1.2, horiz=TRUE, bty='y', inset=c(0,0))
dev.off()


########################################################
#                    MAKE TABLE 5
########################################################
# Test various models
y_m_id <- list() # ss304, nickel and uranium
for(jj in seq_along(m_id)){
  j <- m_id[jj]
  tmp <- read.csv(paste0("data/", fnames[j]))
  X <- tmp[, 12:17]
  ord <- order(X[,1])
  X <- X[ord,]
  for(i in 1:ncol(X)) X[,i] <- BASS:::scale.range(X[,i])
  y_m_id[[jj]] <- tmp[ord,11]
}

X_u <- as.matrix(X)%*%matrix(eigen(C_u)$vectors[,1:2], nrow=6)
mod_us <- bass(X_u, y_m_id[[1]])
#mod_un <- bass(X_u, y_m_id[[2]])
mod_uu <- bass(X_u, y_m_id[[3]])

X_s <- as.matrix(X)%*%matrix(eigen(C_s)$vectors[,1:2], nrow=6)
mod_ss <- bass(X_s, y_m_id[[1]])
#mod_sn <- bass(X_s, y_m_id[[2]])
mod_su <- bass(X_s, y_m_id[[3]])

X_su <- as.matrix(X)%*%matrix(eigen(V_su)$vectors[,c(1, 6)], nrow=6)
mod_sus <- bass(X_su, y_m_id[[1]])
#mod_sn <- bass(X_su, y_m_id[[2]])
mod_suu <- bass(X_su, y_m_id[[3]])

X_zham <- as.matrix(X)%*%matrix(eigen(C_u+C_s)$vectors[,1:2], nrow=6)
mod_z <- bass(X_zahm, y_m_id[[1]])
#mod_sn <- bass(X_su, y_m_id[[2]])
mod_z <- bass(X_zahm, y_m_id[[3]])
