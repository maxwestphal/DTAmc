#-------------------------------------------------------------------------------------------
# Master Thesis: "Adjusting for multiplicity in diagnostic studies: A Wild Bootstrap
# approach for obtaining simultaneous confidence intervals for sensitivity and specificity"
# ------------------------------------------------------------------------------------------
# Author: Anja Rudolph
# Supervisor: Antonia Zapf
# MSc Medical Biometry, Heidelberg
# ------------------------------------------------------------------------------------------
# Program to calculate unadjusted and multiplicity adjusted confidence intervals
# for sensitivity and specificity
# Version 1.0 - last edited 26/02/2017
# Load required libraries
library(multcomp)
library(OptimalCutpoints)
# Function for logit-transformation
logit <- function(p){
  log(p/(1-p))
}
# Function to calculate transformed and backtransformed CIs
ki <- function(est,quant,V){
  auc_logit=log(est/(1-est)) # logit transformation of AUC/SP/SE
  term_l=auc_logit-quant*sqrt(diag(V)) # calculate transformed lower CI
  ki_l=exp(term_l)/(1+exp(term_l)) # backtransformation of lower CI
  return(ki_l)
}
##### ADDED
data <- readr::read_csv2(normalizePath("C:\\Users\\maxwe\\Documents\\R\\PUBLIC\\131022 Daten Ballmann reduziert.csv", "\\"))
data
?normalizePath
d <- 5
n0 <- 80
n1 <- 20
N <- n0 + n1
status <- NA
tag.healthy <- 0
#data <- NA
datatype <- "Example"
datatype <- "Simulated"
alpha <- 0.05
minAUC <- NA
minSE <- 0.9
minSP <- 0.8

?optimal.cutpoints
?control.cutpoints

##### ADDED END
Calc_CI <- function(d,n0,n1,N,status,tag.healthy,data,datatype,alpha,minAUC,minSE,minSP){
  # Results memory----------------------------------------------------------#
  if(datatype == "Example") {
    erg=list()
  }
  if(datatype == "Simulated") {
    ressim=list()
  }
  #-------------------------------------------------------------------------#
  cut <- rep(0,d)
  if(datatype == "Example") {
    for(i in 1:d){
      # Define cut-off depending on pre-specified SE and SP values
      coff <- optimal.cutpoints(X=paste("V",i,sep=""),status,tag.healthy,methods="MinValueSpSe",
                                data,control=control.cutpoints(valueSe=minSE,valueSp=minSP))
      cut[i] <- coff$MinValueSpSe$Global$optimal.cutoff$cutoff
      # To compute sensitivity (SE), replace values of healthy subjects with cut-off values
      data[1:n0,(d+i)] <- cut[i]
      # To compute specificity (SP), replace values of diseased subjects with cut-off values
      data[(n0+1):N,(2*d+i)] <- cut[i]
    }
  }
  # Define partial data sets to compute internal ranks
  x0 <- data[1:n0,1:(3*d)]
  x1 <- data[(n0+1):N,1:(3*d)]
  x <- data[,1:(3*d)]
  # Assign ranks
  rx0 <- apply(x0,2,rank) # internal ranks in healthy
  rx1 <- apply(x1,2,rank) # internal ranks in diseased
  rx <- apply(x,2,rank) # ranks overall
  # Point estimate of AUC/SP/SE computed using matrices
  Pn0 <- diag(n0) - 1/n0
  Pn1 <- diag(n1) - 1/n1
  pl0 <- 1/n1*(rx[1:n0,]-rx0) # also called "normed placements" Z0 in Zapf et al. 2015
  pl1 <- 1/n0*(rx[(n0+1):N,]-rx1) # also called "normed placements" Z1 in Zapf et al. 2015
  l0Z <- Pn0%*%pl0 # multiply Z0 with matrix Pn0
  pl1Z <- Pn1%*%pl1 # multiply Z1 with matrix Pn1
  # pd = point estimate AUC/SP/SE
  pd.orig <- pd <- colMeans(pl1)
  if(datatype == "Example") {
    # Write to memory
    erg[["p_estimate"]] <- pd
  }
  if(datatype == "Simulated") {
    # Calculate bias of point estimators
    bias.se <- (pd[(d+1):(2*d)]-minSE)
    # Write to memory
    ressim[["p_estimate"]] <- pd[(d+1):(2*d)]
    ressim[["bias"]] <- bias.se
  }
  pd0 <- (pd==0) # logical indictor whether AUC/SP/SE = 0
  pd1 <- (pd==1) # logical indictor whether AUC/SP/SE = 1
  # Check for AUC == 0 or 1
  for(i in 1:d){
    if(pd1[[i]] == 1 || pd0[[i]] == 1) {
      if(datatype == "Example") {
        # Replace largest value of healthy subjects with smallest value of diseased
        # subjects to avoid complete segregation of cohorts
        pos=which(data[1:n0,i]==max(data[1:n0,i]))
        data[pos,i]=min(data[(n0+1):N,i])
      }
      if(datatype == "Simulated") {
        # Replace the first largest value of healthy subjects with 1
        pos=which.max(data[1:n0,i])
        data[pos,i]=1
      }
    }
  }
  # Check for SE == 0 or 1
  for(i in (d+1):(2*d)){
    if(pd1[[i]] == 1 || pd0[[i]] == 1) {
      if(datatype == "Example") {
        # Replace largest value of diseased subjects with smallest value of healthy
        # subjects to avoid complete segregation of cohorts
        pos=which(data[(n0+1):N,i]==max(data[(n0+1):N,i]))
        data[pos,i]=min(data[1:n0,i])
      }
      if(datatype == "Simulated") {
        # Replace the first largest value of diseased subjects with 0.5
        pos=which.max(data[(n0+1):N,i])
        data[(n0+pos),i]=0.5
      }
    }
  }
  # Check for SP == 0 or 1
  for(i in (2*d+1):(3*d)){
    if(pd1[[i]] == 1 || pd0[[i]] == 1) {
      if(datatype == "Example") {
        # Replace largest value of healthy subjects with smallest value of diseased
        # subjects to avoid complete segregation of cohorts
        pos=which(data[1:n0,i]==max(data[1:n0,i]))
        data[pos,i]=min(data[(n0+1):N,i])
      }
      if(datatype == "Simulated") {
        # Replace the first largest value of healthy subjects with 0.5
        pos=which.max(data[1:n0,i])
        data[pos,i]=0.5
      }
    }
  }
  # Re-define partial data sets to compute variances
  x0 <- data[1:n0,1:(3*d)]
  x1 <- data[(n0+1):N,1:(3*d)]
  x <- data[,1:(3*d)]
  # Assign ranks
  rx0 <- apply(x0,2,rank) # internal ranks in healthy
  rx1 <- apply(x1,2,rank) # internal ranks in diseased
  rx <- apply(x,2,rank) # ranks overall
  # Point estimate of AUC/SE/SP computed using matrices
  Pn0 <- diag(n0) - 1/n0
  Pn1 <- diag(n1) - 1/n1
  pl0 <- 1/n1*(rx[1:n0,]-rx0) # also called "normed placements" Z0 in Zapf et al. 2015
  pl1 <- 1/n0*(rx[(n0+1):N,]-rx1) # also called "normed placements" Z1 in Zapf et al. 2015
  pl0Z <- Pn0%*%pl0 # multiply Z0 with matrix Pn0
  pl1Z <- Pn1%*%pl1 # multiply Z1 with matrix Pn1
  # pd = point estimate AUC/SE/SP
  pd <- colMeans(pl1)
  # Variance estimators - covariance matrices - multiplying by N is left out here, simplifying
  # calculation of test statistics and confidence intervals further downstream
  VD.auc <- var(pl0[,1:d])/n0 + var(pl1[,1:d])/n1
  VD.se <- var(pl0[,(d+1):(2*d)])/n0 + var(pl1[,(d+1):(2*d)])/n1
  VD.sp <- var(pl0[,(2*d+1):(3*d)])/n0 + var(pl1[,(2*d+1):(3*d)])/n1
  # Scale covariance matrix into correlation matrix
  R.auc <- cov2cor(VD.auc)
  R.se <- cov2cor(VD.se)
  R.sp <- cov2cor(VD.sp)
  # Test statistics
  TStat.auc <- (pd[1:d]-minAUC)/sqrt(c(diag(VD.auc)))
  TStat.se <- (pd[(d+1):(2*d)]-minSE)/sqrt(c(diag(VD.se)))
  TStat.sp <- (pd[(2*d+1):(3*d)]-minSP)/sqrt(c(diag(VD.sp)))
  # Maximal test statistics
  Tmax.auc <- max(TStat.auc)
  Tmax.se <- max(TStat.se)
  Tmax.sp <- max(TStat.sp)
  # Write to memory
  if(datatype == "Example") {
    erg[["TStat"]] <- c(TStat.auc,TStat.se,TStat.sp)
  }
  if(datatype == "Simulated") {
    ressim[["TStat"]] <- TStat.se
  }
  # Logit-transformation
  pdlogit <- logit(pd)
  # SN = Covariance matrices multiplied with Jacobian matrices
  VDlogit.auc <- diag(1/(pd[1:d] *(1-pd[1:d])))%*% VD.auc%*%diag(1/(pd[1:d]* (1-pd[1:d])))
  VDlogit.se <- diag(1/(pd[(d+1):(2*d)] *(1-pd[(d+1):(2*d)])))%*% VD.se%*% diag(1/(pd[(d+1):(2*d)]* (1-pd[(d+1):(2*d)])))
  VDlogit.sp <- diag(1/(pd[(2*d+1):(3*d)]*(1-pd[(2*d+1):(3*d)])))%*%VD.sp%*% diag(1/(pd[(2*d+1):(3*d)]*(1-pd[(2*d+1):(3*d)])))
  # Logit transformed test statistics
  Tlogit.auc <- (pdlogit[1:d]-logit(minAUC))/sqrt(c(diag(VDlogit.auc)))
  Tlogit.se <- (pdlogit[(d+1):(2*d)]-logit(minSE))/sqrt(c(diag(VDlogit.se)))
  Tlogit.sp <- (pdlogit[(2*d+1):(3*d)]-logit(minSP))/sqrt(c(diag(VDlogit.sp)))
  # Maximal test statistics
  Tmlogit.auc <- max(Tlogit.auc)
  Tmlogit.se <- max(Tlogit.se)
  Tmlogit.sp <- max(Tlogit.sp)
  # Write to memory
  if(datatype == "Example") {
    erg[["TLogit"]] <- c(Tlogit.auc,Tlogit.se,Tlogit.sp)
  }
  if(datatype == "Simulated") {
    ressim[["TLogit"]] <- Tlogit.se
  }
  # Define critical values MCP
  qMCP.auc=qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmax.auc,corr=R.auc)$quantile
  qMCP.se =qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmax.se,corr=R.se)$quantile
  qMCP.sp =qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmax.sp,corr=R.sp)$quantile
  if(datatype == "Simulated") {
    # Check whether TStat.se is larger than critical value
    MCP.se <- (TStat.se > qMCP.se)
    MCP <- (mean(MCP.se == 1)) != 0
    # Write to memory
    ressim[["MCP"]] <- c(MCP)
  }
  # Define critical values Logit
  qlog.auc=qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmlogit.auc,corr=R.auc)$quantile
  qlog.se =qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmlogit.se,corr=R.se)$quantile
  qlog.sp =qmvnorm(p=1-(alpha/2),tail="lower",upper=Tmlogit.sp,corr=R.sp)$quantile
  if(datatype == "Simulated") {
    # Check whether Tlogit.se is larger than critical value
    Logit.se <- (Tlogit.se > qlog.se)
    Logit <- (mean(Logit.se == 1)) != 0
    # Write to memory
    ressim[["Logit"]] <- c(Logit)
  }
  # Define critical values unadjusted & Bonferroni
  crit <- qnorm((1-(alpha/2)))
  critBonf <- qnorm(1-(alpha/(2*d)))
  # Unadjusted -----------------------------------------------------------#
  if(datatype == "Simulated") {
    # Check whether TStat.se is larger than critical value
    Unadj.se <- (TStat.se > crit)
    Unadj <- (mean(Unadj.se == 1)) != 0
    # Write to memory
    ressim[["Unadj"]] <- c(Unadj)
    # Unadjusted logit-transformed------------------------------------------#
    # Check whether Tlogit.se is larger than critical value
    Unadjt.se <- (Tlogit.se > crit)
    Unadjt <- (mean(Unadjt.se == 1)) != 0
    # Write to memory
    ressim[["Unadjt"]] <- c(Unadjt)
    # Bonferroni-adjusted --------------------------------------------------#
    # Check whether TStat.se is larger than critical value
    Bonf.se <- (TStat.se > critBonf)
    Bonf <- (mean(Bonf.se == 1)) != 0
    # Write to memory
    ressim[["Bonf"]] <- Bonf
    # Bonferroni-adjusted logit-transformed---------------------------------#
    # Check whether Tlogit.se is larger than critical value
    Bonft.se <- (Tlogit.se > critBonf)
    Bonft <- (mean(Bonft.se == 1)) != 0
    # Write to memory
    ressim[["Bonft"]] <- Bonft
  }
  if(datatype == "Example") {
    # Confidence intervals for comparison (MCP, Logit, unadjusted and Bonferroni)
    # MCP, not logit transformed
    MCP_ki.auc=pd.orig[1:d]-qMCP.auc*sqrt(diag(VD.auc))
    MCP_ki.se =pd.orig[(d+1):(2*d)]-qMCP.se*sqrt(diag(VD.se))
    MCP_ki.sp =pd.orig[(2*d+1):(3*d)]-qMCP.sp*sqrt(diag(VD.sp))
    # Write to memory
    erg[["MCP_ki"]] <- c(MCP_ki.auc,MCP_ki.se,MCP_ki.sp)
    # Logit transformed
    logit_ki.auc=ki(pd[1:d],qlog.auc,VDlogit.auc)
    logit_ki.se =ki(pd[(d+1):(2*d)],qlog.se,VDlogit.se)
    logit_ki.sp =ki(pd[(2*d+1):(3*d)],qlog.sp,VDlogit.sp)
    # Write to memory
    erg[["logit_ki"]] <- c(logit_ki.auc,logit_ki.se,logit_ki.sp)
    # Unadjusted, not logit transformed
    oa_ki.auc=pd.orig[1:d]-qnorm(1-(alpha/2))*sqrt(diag(VD.auc))
    oa_ki.se =pd.orig[(d+1):(2*d)]-qnorm(1-(alpha/2))*sqrt(diag(VD.se))
    oa_ki.sp =pd.orig[(2*d+1):(3*d)]-qnorm(1-(alpha/2))*sqrt(diag(VD.sp))
    # Write to memory
    erg[["oa_ki"]] <- c(oa_ki.auc,oa_ki.se,oa_ki.sp)
    # Unadjusted, logit transformed
    logit_oa_ki.auc=ki(pd[1:d],qnorm(1-(alpha/2)),VDlogit.auc)
    logit_oa_ki.se =ki(pd[(d+1):(2*d)],qnorm(1-(alpha/2)),VDlogit.se)
    logit_oa_ki.sp =ki(pd[(2*d+1):(3*d)],qnorm(1-(alpha/2)),VDlogit.sp)
    # Write to memory
    erg[["logit_oa_ki"]] <- c(logit_oa_ki.auc,logit_oa_ki.se,logit_oa_ki.sp)
    # Bonferroni, not logit transformed
    bonf_ki.auc=pd.orig[1:d]-qnorm(1-(alpha/(2*d)))*sqrt(diag(VD.auc))
    bonf_ki.se =pd.orig[(d+1):(2*d)]-qnorm(1-(alpha/(2*d)))*sqrt(diag(VD.se))
    bonf_ki.sp =pd.orig[(2*d+1):(3*d)]-qnorm(1-(alpha/(2*d)))*sqrt(diag(VD.sp))
    # Write to memory
    erg[["bonf_ki"]] <- c(bonf_ki.auc,bonf_ki.se,bonf_ki.sp)
    # Bonferroni, logit transformed
    logit_bonf_ki.auc=ki(pd[1:d],qnorm(1-(alpha/(2*d))),VDlogit.auc)
    logit_bonf_ki.se =ki(pd[(d+1):(2*d)],qnorm(1-(alpha/(2*d))),VDlogit.se)
    logit_bonf_ki.sp =ki(pd[(2*d+1):(3*d)],qnorm(1-(alpha/(2*d))),VDlogit.sp)
    # Write to memory
    erg[["logit_bonf_ki"]] <- c(logit_bonf_ki.auc,logit_bonf_ki.se,logit_bonf_ki.sp)
  }
  # Wild-Bootstrap approach
  # Create large matrices to be filled
  WR0 = WN0 = WU0 = matrix(0,nrow=nboot, ncol=n0)
  WR1 = WN1 = WU1 = matrix(0,nrow=nboot, ncol=n1)
  for(hh in 1:nboot){
    # Sample from uniform distribution --> Rademacher weights
    WR0[hh, ] <- runif(n0)
    # Sample from normal distribution N(0,1)
    WN0[hh, ] <- rnorm(n0)
    # Sample from uniform distribution between -sqrt(12)/2 and sqrt(12)/2
    WU0[hh, ] <- runif(n0,-sqrt(12)/2,sqrt(12)/2)
    WR1[hh, ] <- runif(n1)
    WN1[hh, ] <- rnorm(n1)
    WU1[hh, ] <- runif(n1,-sqrt(12)/2,sqrt(12)/2)
  }
  WR00<-(WR0<1/2)
  WR01<-(WR0>1/2)
  WR0[WR00] <- -1 # Create Rademacher weights for healthy
  WR0[WR01] <- 1
  WR10 <- (WR1<1/2)
  WR11 <- (WR1>1/2)
  WR1[WR10] <- -1 # Create Rademacher weights for diseased
  WR1[WR11] <- 1
  # Matrices for means
  WRM0 <- -1/n0*WR0
  WRM1 <- 1/n1*WR1
  WNM0 <- -1/n0*WN0
  WNM1 <- 1/n1*WN1
  WUM0 <- -1/n0*WU0
  WUM1 <- 1/n1*WU1
  # Matrices for variances
  WR02 <- WR0^2
  WR12 <- WR1^2
  WN02 <- WN0^2
  WN12 <- WN1^2
  WU02 <- WU0^2
  WU12 <- WU1^2
  WRV0 <- 1/(n0-1)*WR02
  WRV1 <- 1/(n1-1)*WR12
  WNV0 <- 1/(n0-1)*WN02
  WNV1 <- 1/(n1-1)*WN12
  WUV0 <- 1/(n0-1)*WU02
  WUV1 <- 1/(n1-1)*WU12
  # Rademacher weights -----------------------------------------------------------#
  # Multiply Z0 with R. weights
  MZR0 <- WRM0%*%pl0Z
  # Multiply Z1 with R. weights
  MZR1 <- WRM1%*%pl1Z
  # Empirical variance for Z0s
  VZR0 <- (WRV0%*%pl0Z^2-n0*MZR0^2/(n0-1))/(n0)
  # Empirical variance for Z1s
  VZR1 <- (WRV1%*%pl1Z^2-n1*MZR1^2/(n1-1))/(n1)
  # Resampling distribution
  TR <- (MZR0 - MZR1) / sqrt(VZR0 + VZR1)
  QTR.auc=quantile(TR[,1:d],prob=(1-(alpha/2)))
  QTR.se =quantile(TR[,(d+1):(2*d)],prob=(1-(alpha/2)))
  QTR.sp =quantile(TR[,(2*d+1):(3*d)],prob=(1-(alpha/2)))
  if(datatype == "Simulated") {
    # Check whether TStat.se is larger than critical value
    WRade.se <- (TStat.se > QTR.se)
    WRade <- (mean(WRade.se == 1)) != 0
    # Write to memory
    ressim[["WRade"]] <- c(WRade)
  }
  # N(0,1)-weights ---------------------------------------------------------------#
  MZN0 <- WNM0%*%pl0Z
  MZN1 <- WNM1%*%pl1Z
  VZN0 <- (WNV0%*%pl0Z^2-n0*MZN0^2/(n0-1))/(n0)
  VZN1 <- (WNV1%*%pl1Z^2-n1*MZN1^2/(n1-1))/(n1)
  TN <- (MZN0 + MZN1) / sqrt(VZN0 + VZN1)
  QTN.auc=quantile(TN[,1:d],prob=(1-(alpha/2)))
  QTN.se =quantile(TN[,(d+1):(2*d)],prob=(1-(alpha/2)))
  QTN.sp =quantile(TN[,(2*d+1):(3*d)],prob=(1-(alpha/2)))
  if(datatype == "Simulated") {
    # Check whether TStat.se is larger than critical value
    WNormal.se <- (TStat.se > QTN.se)
    WNormal <- (mean(WNormal.se == 1)) != 0
    # Write to memory
    ressim[["WNormal"]] <- c(WNormal)
  }
  # Uniform weights --------------------------------------------------------------#
  MZU0 <- WUM0%*%pl0Z
  MZU1 <- WUM1%*%pl1Z
  VZU0 <- (WUV0%*%pl0Z^2-n0*MZU0^2/(n0-1))/(n0)
  VZU1 <- (WUV1%*%pl1Z^2-n1*MZU1^2/(n1-1))/(n1)
  TU <- (MZU0 + MZU1) / sqrt(VZU0 + VZU1)
  QTU.auc=quantile(TU[,1:d],prob=(1-(alpha/2)))
  QTU.se =quantile(TU[,(d+1):(2*d)],prob=(1-(alpha/2)))
  QTU.sp =quantile(TU[,(2*d+1):(3*d)],prob=(1-(alpha/2)))
  if(datatype == "Simulated") {
    # Check whether TStat.se is larger than critical value
    WUniform.se <- (TStat.se > QTU.se)
    WUniform <- (mean(WUniform.se == 1)) != 0
    # Write to memory
    ressim[["WUniform"]] <- c(WUniform)
  }
  if(datatype == "Example") {
    # Calculate CI based on Wild bootstrap -----------------------------------------#
    tr_ki.auc=pd.orig[1:d]-QTR.auc*sqrt(diag(VD.auc))
    tr_ki.se =pd.orig[(d+1):(2*d)]-QTR.se*sqrt(diag(VD.se))
    tr_ki.sp =pd.orig[(2*d+1):(3*d)]-QTR.sp*sqrt(diag(VD.sp))
    # Write to memory
    erg[["tr_ki"]] <- c(tr_ki.auc,tr_ki.se,tr_ki.sp)
    tn_ki.auc=pd.orig[1:d]-QTN.auc*sqrt(diag(VD.auc))
    tn_ki.se =pd.orig[(d+1):(2*d)]-QTN.se*sqrt(diag(VD.se))
    tn_ki.sp =pd.orig[(2*d+1):(3*d)]-QTN.sp*sqrt(diag(VD.sp))
    # Write to memory
    erg[["tn_ki"]] <- c(tn_ki.auc,tn_ki.se,tn_ki.sp)
    tu_ki.auc=pd.orig[1:d]-QTU.auc*sqrt(diag(VD.auc))
    tu_ki.se =pd.orig[(d+1):(2*d)]-QTU.se*sqrt(diag(VD.se))
    tu_ki.sp =pd.orig[(2*d+1):(3*d)]-QTU.sp*sqrt(diag(VD.sp))
    # Write to memory
    erg[["tu_ki"]] <- c(tu_ki.auc,tu_ki.se,tu_ki.sp)
  }
  # ------------------------------------------------------------------------------#
  if(datatype == "Example") {
    return(erg)
  }
  if(datatype == "Simulated") {
    return(ressim)
  }
}
