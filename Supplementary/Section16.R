############################################
# Last update : October 02 2021
# Code for the simulation in Section 1.6 in the supplementary material
############################################

############################################
# Packages & Source File
############################################

library(lme4)
source("Function16.R")

############################################
# Code
# Note: Recommend to use parallel computing 
# by setting a smaller number of iter.max
############################################

iter.max <- 2                       # Maximum Iteration ; recommended to use a smaller number
Result <- matrix(0,iter.max,66)     # Result matrix
# TTT <- matrix(0,iter.max,4)
# 1:6 columns contain bias estimates of DE under each model specification scenario
# 7:12 columns contain bias estimates of IE under each model specification scenario
# 13:18 columns contain SE estimates of DE under each model specification scenario
# 19:24 columns contain SE estimates of IE under each model specification scenario

for(iteration in 1:iter.max){
    
    ############################
    # Parameters and DGP
    ############################
    
    # set.seed(iteration)
    
    N <- 2000                      # Number of Clusters
    K <- 2                          # Number of Cluster Types
    Mk <- c(3,4)                    # Number of Units in Each Cluster
    pk <- c(0.75,0.25)              # Cluster Type Probability
    
    SigmaX <- 1                     # Variance of Covariates, X ~ N(0,SigmaX^2)
    rhoX <- 0                       # Covariance of Covariates, Cov(X_{ij},X_{ij'}) = rhoX
    
    betaXg <- list()                              
    
    # beta_{g,1} : Outcome Regression Coefficient of Cluster Type 1
    betaXg[[1]] <- c(2,                         # intercept
                     3, 0.8,                    # A_{ij}, A_{i(-j)}
                     1, 0.5,                    # interaction: A_{ij}*C_{i} , A_{i(-j)}*C_{i}
                     -1,0.5,-0.3,0.15,
                     0.8)        # covariate X_{ij}, X_{i(-j)}
    
    # beta_{g,2} : Outcome Regression Coefficient of Cluster Type 2
    betaXg[[2]] <- c(1,
                     2, 0.4,
                     0.5, 0.3,
                     -0.8,0.4,-0.2,0.1,
                     0.6)           
    
    rho <- c(0.1,0.1)               # Variance of Random Effects
    eta <- c(1,1)                   # Variance of Individual Level Error ($\eta$ in the paper is the reciprocal value of eta used in the)
    
    betaXe <- list()
    betaXe[[1]] <- c(-1.25,
                     2, 0.3,
                     0.2, 0.1)     # beta_{e,1} : Propensity Score Coefficient of Cluster Type 1
    betaXe[[2]] <- c(-1,
                     1.25, 0.2,
                     0.15, 0.1)          # beta_{e,2} : Propensity Score Coefficient of Cluster Type 2
    
    lambda <- c(0.25,0.25)                # Variance of Random Effects ($\lambda$ in the paper is the reciprocal value of lambda used in the code)
    
    alpha <- 0.4                    # Policy Parameter of DE
    alpha1 <- 0.8                   # First Policy Parameter of IE
    alpha2 <- 0.2                   # Second Policy Parameter of IE
    
    TrueEffect <- list()                                # True effect list
    TrueEffect[[1]] <- FindEffect(betaXg[[1]],Mk[1],alpha,alpha1,alpha2)    # True DE/IE of cluster type 1
    TrueEffect[[2]] <- FindEffect(betaXg[[2]],Mk[2],alpha,alpha1,alpha2)    # True DE/IE of cluster type 2
    TrueDE <- TrueEffect[[1]]$DE*pk[1] + TrueEffect[[2]]$DE*pk[2]           # True DE
    TrueIE <- TrueEffect[[1]]$IE*pk[1] + TrueEffect[[2]]$IE*pk[2]           # True IE
    
    L <- sort(sample(K,N,replace=TRUE,prob=pk))         # Cluster Type Variable
    Nk <- c(sum(L==1),sum(L==2))                        # Number of Clusters under Type k
    
    Cluster <- list()                                   # Cluster ID Enumeration
    Cluster[[1]] <- rep((1:Nk[1]),each=Mk[1])           # Cluster type 1 ID Enumeration
    Cluster[[2]] <- rep((1:Nk[2]),each=Mk[2])           # Cluster type 2 ID Enumeration
    Cluster[[11]] <- rep(1:N,c(rep(Mk[1],Nk[1]),rep(Mk[2],Nk[2])))       # Cluster type ID Enumeration
    
    X <- list()                                                         # Pre-treatment Covariate List
    X[[1]] <- matrix(0,Nk[1],Mk[1]*3)                                   # N_1 x 4M_1 matrix
    X[[1]][,(1:Mk[1])*3-2] <- matrix(rbinom(Nk[1]*Mk[1],1,0.5),Nk[1],Mk[1])
    X[[1]][,(1:Mk[1])*3  ] <- matrix(rep(rnorm(Nk[1]),each=Mk[1]),Nk[1],Mk[1],byrow=T)
    X[[1]][,(1:Mk[1])*3-1] <- GenerateX(Nk[1],Mk[1],SigmaX,rhoX,X[[1]][,3])
    
    
    X[[2]] <- matrix(0,Nk[2],Mk[2]*3)                                   # N_2 x 4M_2 matrix
    X[[2]][,(1:Mk[2])*3-2] <- matrix(rbinom(Nk[2]*Mk[2],1,0.5),Nk[2],Mk[2])
    X[[2]][,(1:Mk[2])*3  ] <- matrix(rep(rnorm(Nk[2]),each=Mk[2]),Nk[2],Mk[2],byrow=T)
    X[[2]][,(1:Mk[2])*3-1] <- GenerateX(Nk[2],Mk[2],SigmaX,rhoX,X[[2]][,3])
    
    
    WrongX <- list()
    WrongX[[1]] <- matrix(0,Nk[1],Mk[1]*3)
    WrongX[[1]][,(1:Mk[1])*3-2] <- X[[1]][,(1:Mk[1])*3-2]*exp(X[[1]][,(1:Mk[1])*3-1]/4)
    WrongX[[1]][,(1:Mk[1])*3-1] <- X[[1]][,(1:Mk[1])*3-1]^2
    WrongX[[1]][,(1:Mk[1])*3  ] <- X[[1]][,(1:Mk[1])*3-2]*X[[1]][,(1:Mk[1])*3]
    
    WrongX[[2]] <- matrix(0,Nk[2],Mk[2]*3)
    WrongX[[2]][,(1:Mk[2])*3-2] <- X[[2]][,(1:Mk[2])*3-2]*exp(X[[2]][,(1:Mk[2])*3-1]/4)
    WrongX[[2]][,(1:Mk[2])*3-1] <- X[[2]][,(1:Mk[2])*3-1]^2
    WrongX[[2]][,(1:Mk[2])*3  ] <- X[[2]][,(1:Mk[2])*3-2]*X[[2]][,(1:Mk[2])*3]
    
    
    LongX <- list()                                     # Pre-treatment Covariate List (Long Format)
    LongX[[1]] <- ReformX(X[[1]],2,1)                     # N_1*M_1 x 2 matrix
    LongX[[2]] <- ReformX(X[[2]],2,1)                     # N_2*M_2 x 2 matrix
    LongX[[11]] <- rbind(LongX[[1]],LongX[[2]])         # (N_1*M_1+N_2*M_2) x 2 matrix
    
    WrongLongX <- list()
    WrongLongX[[1]] <- ReformX(WrongX[[1]],2,1)                     # N_1*M_1 x 2 matrix
    WrongLongX[[2]] <- ReformX(WrongX[[2]],2,1)                     # N_2*M_2 x 2 matrix
    
    
    position.of.interaction <- c(5,10)
    position.of.X.OR <- 1:5
    position.of.X <- 1:4
    
    A <- list() ; LongA <- list()                       # A list / A list  (Long Format)
    Atemp1 <- GenerateA(LongX[[1]][,position.of.X],Mk[1],betaXe[[1]],lambda[1])
    Atemp2 <- GenerateA(LongX[[2]][,position.of.X],Mk[2],betaXe[[2]],lambda[2])
    A[[1]] <- Atemp1$AShort                             # N_1 x M_1 matrix
    A[[2]] <- Atemp2$AShort                             # N_2 x M_2 matrix
    LongA[[1]] <- Atemp1$ALong                          # N_1*M_1 x 1 matrix
    LongA[[2]] <- Atemp2$ALong                          # N_2*M_2 x 1 matrix
    LongA[[11]] <- rbind(LongA[[1]],LongA[[2]])         # (N_1*M_1+N_2*M_2) x 1 matrix
    
    LongAX <- list()
    # (X_{ij},X_(-{ij}),C_i)*A_{ij} + (X_{ij},X_(-{ij}),C_i)*A_{i(-j)}
    # (1,2,5,6,9,10,11,12)
    LongAX[[1]] <- cbind( LongX[[1]]*LongA[[1]][,1] , LongX[[1]]*LongA[[1]][,2] )
    LongAX[[2]] <- cbind( LongX[[2]]*LongA[[2]][,1] , LongX[[2]]*LongA[[2]][,2] )
    LongAX[[11]] <- rbind(LongAX[[1]],LongAX[[2]])
    
    WrongLongAX <- list()
    # (X_{ij},X_(-{ij}),C_i)*A_{ij} + (X_{ij},X_(-{ij}),C_i)*A_{i(-j)}
    # (1,2,5,6,9,10,11,12)
    WrongLongAX[[1]] <- cbind( WrongLongX[[1]]*LongA[[1]][,1] , WrongLongX[[1]]*LongA[[1]][,2] )
    WrongLongAX[[2]] <- cbind( WrongLongX[[2]]*LongA[[2]][,1] , WrongLongX[[2]]*LongA[[2]][,2] )
    WrongLongAX[[11]] <- rbind(WrongLongAX[[1]],WrongLongAX[[2]])
    
    Y <- list(); LongY <- list()                        # Y list / Y list (Long Format)
    G <- list(); LongG <- list()                        # E(Y|a,x,k) list (and Long Format) 
    Ytemp1 <- GenerateY(LongA[[1]],LongAX[[1]][,position.of.interaction],LongX[[1]][,position.of.X.OR],Mk[1],betaXg[[1]],rho[1],eta[1])
    Ytemp2 <- GenerateY(LongA[[2]],LongAX[[2]][,position.of.interaction],LongX[[2]][,position.of.X.OR],Mk[2],betaXg[[2]],rho[2],eta[2])
    Y[[1]] <- Ytemp1$YShort                             # N_1 x M_1 matrix
    Y[[2]] <- Ytemp2$YShort                             # N_2 x M_2 matrix
    LongY[[1]] <- Ytemp1$YLong                          # N_1*M_1 x 1 matrix
    LongY[[2]] <- Ytemp2$YLong                          # N_2*M_2 x 1 matrix
    LongY[[11]] <- rbind(LongY[[1]],LongY[[2]])         # (N_1*M_1+N_2*M_2) x 1 matrix
    
    G[[1]] <- Ytemp1$YmeanShort                         # N_1 x M_1 matrix
    G[[2]] <- Ytemp2$YmeanShort                         # N_2 x M_2 matrix
    LongG[[1]] <- Ytemp1$YmeanLong                      # N_1*M_1 x 1 matrix
    LongG[[2]] <- Ytemp2$YmeanLong                      # N_2*M_2 x 1 matrix
    LongG[[11]] <- rbind(LongG[[1]],LongG[[2]])         # (N_1*M_1+N_2*M_2) x 1 matrix
    
    ############################
    # Over-specified Groups
    ############################
    
    wK <- 4                         # Number of Over-specified Cluster Types
    wL <- rep(0,N)                  # f-specified Cluster Type Variable
    for(ii in 1:Nk[1]){
        wL[ii] <- ifelse( X[[1]][ii, 3]<1.5 , 1, 2 )
    }
    for(ii in 1:Nk[2]){
        wL[ii+Nk[1]] <- ifelse( X[[2]][ii, 3] < 1.5, 3, 4 )
    }
    wGP <- rep(wL,rep(Mk,Nk))       # Wrong GP
    wL.position <- list()
    wL.position[[1]] <- which(wL==1)
    wL.position[[2]] <- which(wL==2)
    wL.position[[3]] <- which(wL==3)-Nk[1]
    wL.position[[4]] <- which(wL==4)-Nk[1]
    
    wNk <- c(sum(wL==1),sum(wL==2),sum(wL==3),sum(wL==4))
    
    ############################
    # Estimation under 
    # four model specification
    ############################
    
    # under.position.of.interaction <- c(5,10)
    under.position.of.interaction <- NULL
    under.position.of.X <- c(1:4)
    
    over.position.of.interaction <- 1:10
    over.position.of.X <- c(1:5)
    
    EstG <- list()                                      # Outcome Regression Estimation
    EstG[[1]] <- lmer(LongY[[1]][,1]~LongA[[1]]+LongAX[[1]][,position.of.interaction]+
                          LongX[[1]][,position.of.X.OR]+(1|Cluster[[1]]),REML=FALSE)
    EstG[[2]] <- lmer(LongY[[2]][,1]~LongA[[2]]+LongAX[[2]][,position.of.interaction]+
                          LongX[[2]][,position.of.X.OR]+(1|Cluster[[2]]),REML=FALSE)
    
    USEstG <- list()                                 # Mis-specified Outcome Regression Estimation
    USEstG[[1]] <- lmer(LongY[[1]][,1]~LongA[[1]]+ # WrongLongAX[[1]][,under.position.of.interaction]+
                            WrongLongX[[1]][,under.position.of.X]+(1|Cluster[[1]]),REML=FALSE)
    USEstG[[2]] <- lmer(LongY[[2]][,1]~LongA[[2]]+ # WrongLongAX[[2]][,under.position.of.interaction]+
                            WrongLongX[[2]][,under.position.of.X]+(1|Cluster[[2]]),REML=FALSE)
    
    OSEstG <- list()                                 # Mis-specified Outcome Regression Estimation
    OSEstG[[1]] <- lmer(LongY[[1]][,1]~LongA[[1]]+LongAX[[1]][,over.position.of.interaction]+
                            LongX[[1]][,over.position.of.X]+(1|Cluster[[1]]),REML=FALSE)
    OSEstG[[2]] <- lmer(LongY[[2]][,1]~LongA[[2]]+LongAX[[2]][,over.position.of.interaction]+
                            LongX[[2]][,over.position.of.X]+(1|Cluster[[2]]),REML=FALSE)
    
    OverEstG <- list()                                  # Over-specified Outcome Regression Estimaton
    OverEstG[[1]] <- lmer(LongY[[11]][wGP==1,1]~LongA[[11]][wGP==1,]+LongAX[[11]][wGP==1,position.of.interaction]+
                              LongX[[11]][wGP==1,position.of.X]+(1|Cluster[[11]][wGP==1]),REML=FALSE)
    OverEstG[[2]] <- lmer(LongY[[11]][wGP==2,1]~LongA[[11]][wGP==2,]+LongAX[[11]][wGP==2,position.of.interaction]+
                              LongX[[11]][wGP==2,position.of.X]+(1|Cluster[[11]][wGP==2]),REML=FALSE)
    OverEstG[[3]] <- lmer(LongY[[11]][wGP==3,1]~LongA[[11]][wGP==3,]+LongAX[[11]][wGP==3,position.of.interaction]+
                              LongX[[11]][wGP==3,position.of.X]+(1|Cluster[[11]][wGP==3]),REML=FALSE)
    OverEstG[[4]] <- lmer(LongY[[11]][wGP==4,1]~LongA[[11]][wGP==4,]+LongAX[[11]][wGP==4,position.of.interaction]+
                              LongX[[11]][wGP==4,position.of.X]+(1|Cluster[[11]][wGP==4]),REML=FALSE)
    
    UnderEstG <- list()                                  # Under-specified Outcome Regression Estimaton
    UnderEstG[[1]] <- lmer(LongY[[11]][,1]~LongA[[11]]+LongAX[[11]][,position.of.interaction]+LongX[[11]][,position.of.X]+(1|Cluster[[11]]),REML=FALSE)
    UnderEstG[[2]] <- UnderEstG[[1]]
    
    estbetaXg <- list()                                 # beta_g Estimates
    estbetaXg[[1]] <- fixef(EstG[[1]])                  # beta_{g,1} Estimate
    estbetaXg[[2]] <- fixef(EstG[[2]])                  # beta_{g,2} Estimate
    
    estrho <- c( as.numeric(VarCorr(EstG[[1]])),             # rho estimates
                 as.numeric(VarCorr(EstG[[2]])) )
    
    esteta <- c(sigma(EstG[[1]])^2 , sigma(EstG[[2]])^2)     # eta estimates
    
    USestbetaXg <- list()                                 # Mis-specified beta_g Estimates
    USestbetaXg[[1]] <- fixef(USEstG[[1]])             # beta_{g,1} Estimate
    USestbetaXg[[2]] <- fixef(USEstG[[2]])             # beta_{g,2} Estimate
    
    USestrho <- c( as.numeric(VarCorr(USEstG[[1]])),   # Mis-specified rho estimates
                   as.numeric(VarCorr(USEstG[[2]])) )
    
    USesteta <- c(sigma(USEstG[[1]])^2 ,               # Mis-specified eta estimates
                  sigma(USEstG[[2]])^2)
    
    OSestbetaXg <- list()                                 # Mis-specified beta_g Estimates
    OSestbetaXg[[1]] <- fixef(OSEstG[[1]])             # beta_{g,1} Estimate
    OSestbetaXg[[2]] <- fixef(OSEstG[[2]])             # beta_{g,2} Estimate
    
    OSestrho <- c( as.numeric(VarCorr(OSEstG[[1]])),   # Mis-specified rho estimates
                   as.numeric(VarCorr(OSEstG[[2]])) )
    
    OSesteta <- c(sigma(OSEstG[[1]])^2 ,               # Mis-specified eta estimates
                  sigma(OSEstG[[2]])^2)
    
    
    OverestbetaXg <- list()                                 # Over-specified beta_g Estimates
    OverestbetaXg[[1]] <- fixef(OverEstG[[1]])              # beta_{g,1} Estimate
    OverestbetaXg[[2]] <- fixef(OverEstG[[2]])              # beta_{g,2} Estimate
    OverestbetaXg[[3]] <- fixef(OverEstG[[3]])              # beta_{g,3} Estimate
    OverestbetaXg[[4]] <- fixef(OverEstG[[4]])              # beta_{g,4} Estimate
    
    Overestrho <- c( as.numeric(VarCorr(OverEstG[[1]])),    # Over-specified rho estimates
                     as.numeric(VarCorr(OverEstG[[2]])),
                     as.numeric(VarCorr(OverEstG[[3]])),
                     as.numeric(VarCorr(OverEstG[[4]])) )
    
    Overesteta <- c( sigma(OverEstG[[1]])^2 ,               # Over-specified eta estimates
                     sigma(OverEstG[[2]])^2 ,
                     sigma(OverEstG[[3]])^2 ,
                     sigma(OverEstG[[4]])^2 )
    
    
    
    UnderestbetaXg <- list()                                 # Under-specified beta_g Estimates
    UnderestbetaXg[[1]] <- fixef(UnderEstG[[1]])              # beta_{g,1} Estimate
    UnderestbetaXg[[2]] <- fixef(UnderEstG[[1]])              # beta_{g,1} Estimate
    
    Underestrho <- c( as.numeric(VarCorr(UnderEstG[[1]])),    # Under-specified rho estimates
                      as.numeric(VarCorr(UnderEstG[[1]])) )
    
    Underesteta <- c( sigma(UnderEstG[[1]])^2,                # Under-specified eta estimates
                      sigma(UnderEstG[[1]])^2 )
    
    
    EstE <- list()                                           # Propensity Score Estimation
    EstE[[1]] <- glmer(LongA[[1]][,1]~LongX[[1]][,position.of.X]+(1|Cluster[[1]]),family="binomial")
    EstE[[2]] <- glmer(LongA[[2]][,1]~LongX[[2]][,position.of.X]+(1|Cluster[[2]]),family="binomial")
    
    USEstE <- list()                                      # Propensity Score Estimation
    USEstE[[1]] <- glmer(LongA[[1]][,1]~WrongLongX[[1]][,under.position.of.X]+(1|Cluster[[1]]),family="binomial")
    USEstE[[2]] <- glmer(LongA[[2]][,1]~WrongLongX[[2]][,under.position.of.X]+(1|Cluster[[2]]),family="binomial")
    
    OSEstE <- list()                                      # Propensity Score Estimation
    OSEstE[[1]] <- glmer(LongA[[1]][,1]~LongX[[1]][,over.position.of.X]+(1|Cluster[[1]]),family="binomial")
    OSEstE[[2]] <- glmer(LongA[[2]][,1]~LongX[[2]][,over.position.of.X]+(1|Cluster[[2]]),family="binomial")
    
    OverEstE <- list()                                      # Propensity Score Estimation
    OverEstE[[1]] <- glmer(LongA[[11]][wGP==1,1]~LongX[[11]][wGP==1,position.of.X]+(1|Cluster[[11]][wGP==1]),family="binomial")
    OverEstE[[2]] <- glmer(LongA[[11]][wGP==2,1]~LongX[[11]][wGP==2,position.of.X]+(1|Cluster[[11]][wGP==2]),family="binomial")
    OverEstE[[3]] <- glmer(LongA[[11]][wGP==3,1]~LongX[[11]][wGP==3,position.of.X]+(1|Cluster[[11]][wGP==3]),family="binomial")
    OverEstE[[4]] <- glmer(LongA[[11]][wGP==4,1]~LongX[[11]][wGP==4,position.of.X]+(1|Cluster[[11]][wGP==4]),family="binomial")
    
    UnderEstE <- list()                                      # Propensity Score Estimation
    UnderEstE[[1]] <- glmer(LongA[[11]][,1]~LongX[[11]][,position.of.X]+(1|Cluster[[11]]),family="binomial")
    
    estbetaXe <- list()                                 # beta_e Estimates
    estbetaXe[[1]] <- fixef(EstE[[1]])                  # beta_{e,1} Estimate
    estbetaXe[[2]] <- fixef(EstE[[2]])                  # beta_{e,2} Estimate
    estlambda <- c( as.numeric(VarCorr(EstE[[1]])),     # lambda Estimates
                    as.numeric(VarCorr(EstE[[2]])) )
    
    USestbetaXe <- list()                                      # Mis-specified beta_e Estimates
    USestbetaXe[[1]] <- fixef(USEstE[[1]])                  # beta_{e,1} Estimate
    USestbetaXe[[2]] <- fixef(USEstE[[2]])                  # beta_{e,2} Estimate
    USestlambda <- c( as.numeric(VarCorr(USEstE[[1]])),     # Mis-specified lambda Estimates
                      as.numeric(VarCorr(USEstE[[2]])) )
    
    OSestbetaXe <- list()                                      # Mis-specified beta_e Estimates
    OSestbetaXe[[1]] <- fixef(OSEstE[[1]])                  # beta_{e,1} Estimate
    OSestbetaXe[[2]] <- fixef(OSEstE[[2]])                  # beta_{e,2} Estimate
    OSestlambda <- c( as.numeric(VarCorr(OSEstE[[1]])),     # Mis-specified lambda Estimates
                      as.numeric(VarCorr(OSEstE[[2]])) )
    
    OverestbetaXe <- list()                                      # Over-specified beta_e Estimates
    OverestbetaXe[[1]] <- fixef(OverEstE[[1]])                   # beta_{e,1} Estimate
    OverestbetaXe[[2]] <- fixef(OverEstE[[2]])                   # beta_{e,2} Estimate
    OverestbetaXe[[3]] <- fixef(OverEstE[[3]])                   # beta_{e,3} Estimate
    OverestbetaXe[[4]] <- fixef(OverEstE[[4]])                   # beta_{e,4} Estimate
    Overestlambda <- c( as.numeric(VarCorr(OverEstE[[1]])),      # Mis-specified lambda Estimates
                        as.numeric(VarCorr(OverEstE[[2]])),
                        as.numeric(VarCorr(OverEstE[[3]])),
                        as.numeric(VarCorr(OverEstE[[4]])))
    
    UnderestbetaXe <- list()                                      # Over-specified beta_e Estimates
    UnderestbetaXe[[1]] <- fixef(UnderEstE[[1]])                   # beta_{e,1} Estimate
    UnderestbetaXe[[2]] <- fixef(UnderEstE[[1]])
    Underestlambda <- c( as.numeric(VarCorr(UnderEstE[[1]])),      # Mis-specified lambda Estimates
                         as.numeric(VarCorr(UnderEstE[[1]])) )
    
    ############################
    # PS/OR estimates
    ############################
    
    TruePS <- list()                                    # True PS : Integrate w.r.t. Random Effects
    TruePS[[1]] <- rep(0,Nk[1])
    TruePS[[2]] <- rep(0,Nk[2])
    EstPS <- list()                                     # Est PS : Integrate w.r.t. Random Effects
    EstPS[[1]] <- rep(0,Nk[1])
    EstPS[[2]] <- rep(0,Nk[2])
    USEstPS <- list()                                # Mis-specified Est PS : Integrate w.r.t. Random Effects
    USEstPS[[1]] <- rep(0,Nk[1])
    USEstPS[[2]] <- rep(0,Nk[2])
    OSEstPS <- list()                                # Mis-specified Est PS : Integrate w.r.t. Random Effects
    OSEstPS[[1]] <- rep(0,Nk[1])
    OSEstPS[[2]] <- rep(0,Nk[2])
    OverEstPS <- list()                                # Mis-specified Est PS : Integrate w.r.t. Random Effects
    OverEstPS[[1]] <- rep(0,wNk[1])
    OverEstPS[[2]] <- rep(0,wNk[2])
    OverEstPS[[3]] <- rep(0,wNk[3])
    OverEstPS[[4]] <- rep(0,wNk[4])
    UnderEstPS <- list()                                # Mis-specified Est PS : Integrate w.r.t. Random Effects
    UnderEstPS[[1]] <- rep(0,Nk[1])
    UnderEstPS[[2]] <- rep(0,Nk[2])
    
    
    for(ii in 1:Nk[1]){
        TruePS[[1]][ii] <- truePS(A[[1]][ii,],LongX[[1]][Mk[1]*(ii-1)+(1:Mk[1]),position.of.X],Mk[1],betaXe[[1]],lambda[1])
        EstPS[[1]][ii] <- truePS(A[[1]][ii,],LongX[[1]][Mk[1]*(ii-1)+(1:Mk[1]),position.of.X],Mk[1],estbetaXe[[1]],estlambda[1])
        USEstPS[[1]][ii] <- truePS(A[[1]][ii,],WrongLongX[[1]][Mk[1]*(ii-1)+(1:Mk[1]),under.position.of.X],Mk[1],USestbetaXe[[1]],USestlambda[1])
        OSEstPS[[1]][ii] <- truePS(A[[1]][ii,],LongX[[1]][Mk[1]*(ii-1)+(1:Mk[1]),over.position.of.X],Mk[1],OSestbetaXe[[1]],OSestlambda[1])
        UnderEstPS[[1]][ii] <- truePS(A[[1]][ii,],LongX[[1]][Mk[1]*(ii-1)+(1:Mk[1]),position.of.X],Mk[1],UnderestbetaXe[[1]],Underestlambda[1])
    }
    for(ii in 1:Nk[2]){
        TruePS[[2]][ii] <- truePS(A[[2]][ii,],LongX[[2]][Mk[2]*(ii-1)+(1:Mk[2]),position.of.X],Mk[2],betaXe[[2]],lambda[2])
        EstPS[[2]][ii] <- truePS(A[[2]][ii,],LongX[[2]][Mk[2]*(ii-1)+(1:Mk[2]),position.of.X],Mk[2],estbetaXe[[2]],estlambda[2])
        USEstPS[[2]][ii] <- truePS(A[[2]][ii,],WrongLongX[[2]][Mk[2]*(ii-1)+(1:Mk[2]),under.position.of.X],Mk[2],USestbetaXe[[2]],USestlambda[2])
        OSEstPS[[2]][ii] <- truePS(A[[2]][ii,],LongX[[2]][Mk[2]*(ii-1)+(1:Mk[2]),over.position.of.X],Mk[2],OSestbetaXe[[2]],OSestlambda[2])
        UnderEstPS[[2]][ii] <- truePS(A[[2]][ii,],LongX[[2]][Mk[2]*(ii-1)+(1:Mk[2]),position.of.X],Mk[2],UnderestbetaXe[[2]],Underestlambda[2])
    }
    for(ii in 1:wNk[[1]]){
        OverEstPS[[1]][ii] <- truePS(A[[1]][wL.position[[1]][ii],],LongX[[1]][Mk[1]*(wL.position[[1]][ii]-1)+(1:Mk[1]),position.of.X],Mk[1],OverestbetaXe[[1]],Overestlambda[1])
    }
    for(ii in 1:wNk[[2]]){
        OverEstPS[[2]][ii] <- truePS(A[[1]][wL.position[[2]][ii],],LongX[[1]][Mk[1]*(wL.position[[2]][ii]-1)+(1:Mk[1]),position.of.X],Mk[1],OverestbetaXe[[2]],Overestlambda[2])
    }
    for(ii in 1:wNk[[3]]){
        OverEstPS[[3]][ii] <- truePS(A[[2]][wL.position[[3]][ii],],LongX[[2]][Mk[2]*(wL.position[[3]][ii]-1)+(1:Mk[2]),position.of.X],Mk[2],OverestbetaXe[[3]],Overestlambda[3])
    }
    for(ii in 1:wNk[[4]]){
        OverEstPS[[4]][ii] <- truePS(A[[2]][wL.position[[4]][ii],],LongX[[2]][Mk[2]*(wL.position[[4]][ii]-1)+(1:Mk[2]),position.of.X],Mk[2],OverestbetaXe[[4]],Overestlambda[4])
    }
    
    
    FittedY <- list()                                   # Fitted Y
    FittedY[[1]] <- matrix(cbind(1,LongA[[1]],LongAX[[1]][,position.of.interaction],LongX[[1]][,position.of.X.OR])%*%estbetaXg[[1]],Nk[1],Mk[1],byrow=T)
    FittedY[[2]] <- matrix(cbind(1,LongA[[2]],LongAX[[2]][,position.of.interaction],LongX[[2]][,position.of.X.OR])%*%estbetaXg[[2]],Nk[2],Mk[2],byrow=T)
    
    USFittedY <- list()                              # Mis-specified Fitted Y
    USFittedY[[1]] <- matrix(cbind(1,LongA[[1]],# WrongLongAX[[1]][,under.position.of.interaction],
                                   WrongLongX[[1]][,under.position.of.X])%*%USestbetaXg[[1]],Nk[1],Mk[1],byrow=T)
    USFittedY[[2]] <- matrix(cbind(1,LongA[[2]],# WrongLongAX[[2]][,under.position.of.interaction],
                                   WrongLongX[[2]][,under.position.of.X])%*%USestbetaXg[[2]],Nk[2],Mk[2],byrow=T)
    
    OSFittedY <- list()                              # Mis-specified Fitted Y
    OSFittedY[[1]] <- matrix(cbind(1,LongA[[1]],LongAX[[1]][,over.position.of.interaction],LongX[[1]][,over.position.of.X])%*%OSestbetaXg[[1]],Nk[1],Mk[1],byrow=T)
    OSFittedY[[2]] <- matrix(cbind(1,LongA[[2]],LongAX[[2]][,over.position.of.interaction],LongX[[2]][,over.position.of.X])%*%OSestbetaXg[[2]],Nk[2],Mk[2],byrow=T)
    
    OverFittedY <- list()                              # Over-specified Fitted Y
    OverFittedY[[1]] <- matrix(cbind(1,LongA[[11]][wGP==1,],LongAX[[11]][wGP==1,position.of.interaction],LongX[[11]][wGP==1,position.of.X])%*%OverestbetaXg[[1]],wNk[1],Mk[1],byrow=T)
    OverFittedY[[2]] <- matrix(cbind(1,LongA[[11]][wGP==2,],LongAX[[11]][wGP==2,position.of.interaction],LongX[[11]][wGP==2,position.of.X])%*%OverestbetaXg[[2]],wNk[2],Mk[1],byrow=T)
    OverFittedY[[3]] <- matrix(cbind(1,LongA[[11]][wGP==3,],LongAX[[11]][wGP==3,position.of.interaction],LongX[[11]][wGP==3,position.of.X])%*%OverestbetaXg[[3]],wNk[3],Mk[2],byrow=T)
    OverFittedY[[4]] <- matrix(cbind(1,LongA[[11]][wGP==4,],LongAX[[11]][wGP==4,position.of.interaction],LongX[[11]][wGP==4,position.of.X])%*%OverestbetaXg[[4]],wNk[4],Mk[2],byrow=T)
    
    UnderFittedY <- list()                              # Under-specified Fitted Y
    UnderFittedY[[1]] <- matrix(cbind(1,LongA[[1]],LongAX[[1]][,position.of.interaction],LongX[[1]][,position.of.X])%*%UnderestbetaXg[[1]],Nk[1],Mk[1],byrow=T)
    UnderFittedY[[2]] <- matrix(cbind(1,LongA[[2]],LongAX[[2]][,position.of.interaction],LongX[[2]][,position.of.X])%*%UnderestbetaXg[[2]],Nk[2],Mk[2],byrow=T)
    
    
    Residual <- list()                                  # Residual
    Residual[[1]] <- Y[[1]] - FittedY[[1]]
    Residual[[2]] <- Y[[2]] - FittedY[[2]]
    
    USResidual <- list()                             # Mis-specified Residual
    USResidual[[1]] <- Y[[1]] - USFittedY[[1]]
    USResidual[[2]] <- Y[[2]] - USFittedY[[2]]
    
    OSResidual <- list()                             # Mis-specified Residual
    OSResidual[[1]] <- Y[[1]] - OSFittedY[[1]]
    OSResidual[[2]] <- Y[[2]] - OSFittedY[[2]]
    
    OverResidual <- list()                             # Over-specified Residual
    OverResidual[[1]] <- Y[[1]][wL.position[[1]],] - OverFittedY[[1]]
    OverResidual[[2]] <- Y[[1]][wL.position[[2]],] - OverFittedY[[2]]
    OverResidual[[3]] <- Y[[2]][wL.position[[3]],] - OverFittedY[[3]]
    OverResidual[[4]] <- Y[[2]][wL.position[[4]],] - OverFittedY[[4]]
    
    UnderResidual <- list()                             # Under-specified Residual
    UnderResidual[[1]] <- Y[[1]] - UnderFittedY[[1]]
    UnderResidual[[2]] <- Y[[2]] - UnderFittedY[[2]]
    
    
    ############################
    # Calculate EIF estimates
    ############################
    
    OOCPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    OOCPPart1[[1]] <- ResidualWeighting(OSResidual[[1]],A[[1]],EstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    OOCPPart1[[2]] <- ResidualWeighting(OSResidual[[2]],A[[2]],EstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    COCPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    COCPPart1[[1]] <- ResidualWeighting(Residual[[1]],A[[1]],EstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    COCPPart1[[2]] <- ResidualWeighting(Residual[[2]],A[[2]],EstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    UOCPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    UOCPPart1[[1]] <- ResidualWeighting(USResidual[[1]],A[[1]],EstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    UOCPPart1[[2]] <- ResidualWeighting(USResidual[[2]],A[[2]],EstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    
    
    
    OOOPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    OOOPPart1[[1]] <- ResidualWeighting(OSResidual[[1]],A[[1]],OSEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    OOOPPart1[[2]] <- ResidualWeighting(OSResidual[[2]],A[[2]],OSEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    COOPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    COOPPart1[[1]] <- ResidualWeighting(Residual[[1]],A[[1]],OSEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    COOPPart1[[2]] <- ResidualWeighting(Residual[[2]],A[[2]],OSEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    UOOPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    UOOPPart1[[1]] <- ResidualWeighting(USResidual[[1]],A[[1]],OSEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    UOOPPart1[[2]] <- ResidualWeighting(USResidual[[2]],A[[2]],OSEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    
    
    
    OOUPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    OOUPPart1[[1]] <- ResidualWeighting(OSResidual[[1]],A[[1]],USEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    OOUPPart1[[2]] <- ResidualWeighting(OSResidual[[2]],A[[2]],USEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    COUPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    COUPPart1[[1]] <- ResidualWeighting(Residual[[1]],A[[1]],USEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    COUPPart1[[2]] <- ResidualWeighting(Residual[[2]],A[[2]],USEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    UOUPPart1 <- list()                                   # Residual weighting: [Wrong residual]/[estimated PS]
    UOUPPart1[[1]] <- ResidualWeighting(USResidual[[1]],A[[1]],USEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    UOUPPart1[[2]] <- ResidualWeighting(USResidual[[2]],A[[2]],USEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    OvPart1 <- list()
    OvPart1[[1]] <- ResidualWeighting(OverResidual[[1]],A[[1]][wL.position[[1]],],OverEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    OvPart1[[2]] <- ResidualWeighting(OverResidual[[2]],A[[1]][wL.position[[2]],],OverEstPS[[2]],Mk[1],alpha,alpha1,alpha2)
    OvPart1[[3]] <- ResidualWeighting(OverResidual[[3]],A[[2]][wL.position[[3]],],OverEstPS[[3]],Mk[2],alpha,alpha1,alpha2)
    OvPart1[[4]] <- ResidualWeighting(OverResidual[[4]],A[[2]][wL.position[[4]],],OverEstPS[[4]],Mk[2],alpha,alpha1,alpha2)
    
    UdPart1 <- list()
    UdPart1[[1]] <- ResidualWeighting(UnderResidual[[1]],A[[1]],UnderEstPS[[1]],Mk[1],alpha,alpha1,alpha2)
    UdPart1[[2]] <- ResidualWeighting(UnderResidual[[2]],A[[2]],UnderEstPS[[2]],Mk[2],alpha,alpha1,alpha2)
    
    
    
    
    COPart2 <- list()                                     # Outcome Weighting; sum of {w_k(a)}'\hat{g}(a,x,k)
    COPart2[[1]] <- ORmean(LongX[[1]],position.of.interaction,position.of.X,Mk[1],estbetaXg[[1]],alpha,alpha1,alpha2)
    COPart2[[2]] <- ORmean(LongX[[2]],position.of.interaction,position.of.X,Mk[2],estbetaXg[[2]],alpha,alpha1,alpha2)
    
    UOPart2 <- list()                                     # Outcome Weighting; sum of {w_k(a)}'\hat{g}(a,x,k)
    UOPart2[[1]] <- ORmean(LongX[[1]],under.position.of.interaction,under.position.of.X,Mk[1],USestbetaXg[[1]],alpha,alpha1,alpha2)
    UOPart2[[2]] <- ORmean(LongX[[2]],under.position.of.interaction,under.position.of.X,Mk[2],USestbetaXg[[2]],alpha,alpha1,alpha2)
    
    OOPart2 <- list()                                     # Outcome Weighting; sum of {w_k(a)}'\hat{g}(a,x,k)
    OOPart2[[1]] <- ORmean(LongX[[1]],over.position.of.interaction,over.position.of.X,Mk[1],OSestbetaXg[[1]],alpha,alpha1,alpha2)
    OOPart2[[2]] <- ORmean(LongX[[2]],over.position.of.interaction,over.position.of.X,Mk[2],OSestbetaXg[[2]],alpha,alpha1,alpha2)
    
    OvPart2 <- list()
    OvPart2[[1]] <- ORmean(LongX[[11]][wGP==1,],position.of.interaction,position.of.X,Mk[1],OverestbetaXg[[1]],alpha,alpha1,alpha2)
    OvPart2[[2]] <- ORmean(LongX[[11]][wGP==2,],position.of.interaction,position.of.X,Mk[1],OverestbetaXg[[2]],alpha,alpha1,alpha2)
    OvPart2[[3]] <- ORmean(LongX[[11]][wGP==3,],position.of.interaction,position.of.X,Mk[2],OverestbetaXg[[3]],alpha,alpha1,alpha2)
    OvPart2[[4]] <- ORmean(LongX[[11]][wGP==4,],position.of.interaction,position.of.X,Mk[2],OverestbetaXg[[4]],alpha,alpha1,alpha2)
    
    UdPart2 <- list()
    UdPart2[[1]] <- ORmean(LongX[[1]],position.of.interaction,position.of.X,Mk[1],UnderestbetaXg[[1]],alpha,alpha1,alpha2)
    UdPart2[[2]] <- ORmean(LongX[[2]],position.of.interaction,position.of.X,Mk[2],UnderestbetaXg[[2]],alpha,alpha1,alpha2)
    
    
    ############################
    # Calculate bias/SE
    ############################
    
    COCPDE.cluster <- c( mean(COCPPart1[[1]]$DE + COPart2[[1]]$DE) , mean(COCPPart1[[2]]$DE + COPart2[[2]]$DE) )            # Cluster Level DE
    COCPIE.cluster <- c( mean(COCPPart1[[1]]$IE + COPart2[[1]]$IE) , mean(COCPPart1[[2]]$IE + COPart2[[2]]$IE) )            # Cluster Level IE
    COCPDE.final <- sum(COCPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    COCPIE.final <- sum(COCPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    COOPDE.cluster <- c( mean(COOPPart1[[1]]$DE + COPart2[[1]]$DE) , mean(COOPPart1[[2]]$DE + COPart2[[2]]$DE) )            # Cluster Level DE
    COOPIE.cluster <- c( mean(COOPPart1[[1]]$IE + COPart2[[1]]$IE) , mean(COOPPart1[[2]]$IE + COPart2[[2]]$IE) )            # Cluster Level IE
    COOPDE.final <- sum(COOPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    COOPIE.final <- sum(COOPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    COUPDE.cluster <- c( mean(COUPPart1[[1]]$DE + COPart2[[1]]$DE) , mean(COUPPart1[[2]]$DE + COPart2[[2]]$DE) )            # Cluster Level DE
    COUPIE.cluster <- c( mean(COUPPart1[[1]]$IE + COPart2[[1]]$IE) , mean(COUPPart1[[2]]$IE + COPart2[[2]]$IE) )            # Cluster Level IE
    COUPDE.final <- sum(COUPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    COUPIE.final <- sum(COUPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    
    OOCPDE.cluster <- c( mean(OOCPPart1[[1]]$DE + OOPart2[[1]]$DE) , mean(OOCPPart1[[2]]$DE + OOPart2[[2]]$DE) )            # Cluster Level DE
    OOCPIE.cluster <- c( mean(OOCPPart1[[1]]$IE + OOPart2[[1]]$IE) , mean(OOCPPart1[[2]]$IE + OOPart2[[2]]$IE) )            # Cluster Level IE
    OOCPDE.final <- sum(OOCPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    OOCPIE.final <- sum(OOCPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    OOOPDE.cluster <- c( mean(OOOPPart1[[1]]$DE + OOPart2[[1]]$DE) , mean(OOOPPart1[[2]]$DE + OOPart2[[2]]$DE) )            # Cluster Level DE
    OOOPIE.cluster <- c( mean(OOOPPart1[[1]]$IE + OOPart2[[1]]$IE) , mean(OOOPPart1[[2]]$IE + OOPart2[[2]]$IE) )            # Cluster Level IE
    OOOPDE.final <- sum(OOOPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    OOOPIE.final <- sum(OOOPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    OOUPDE.cluster <- c( mean(OOUPPart1[[1]]$DE + OOPart2[[1]]$DE) , mean(OOUPPart1[[2]]$DE + OOPart2[[2]]$DE) )            # Cluster Level DE
    OOUPIE.cluster <- c( mean(OOUPPart1[[1]]$IE + OOPart2[[1]]$IE) , mean(OOUPPart1[[2]]$IE + OOPart2[[2]]$IE) )            # Cluster Level IE
    OOUPDE.final <- sum(OOUPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    OOUPIE.final <- sum(OOUPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    
    UOCPDE.cluster <- c( mean(UOCPPart1[[1]]$DE + UOPart2[[1]]$DE) , mean(UOCPPart1[[2]]$DE + UOPart2[[2]]$DE) )            # Cluster Level DE
    UOCPIE.cluster <- c( mean(UOCPPart1[[1]]$IE + UOPart2[[1]]$IE) , mean(UOCPPart1[[2]]$IE + UOPart2[[2]]$IE) )            # Cluster Level IE
    UOCPDE.final <- sum(UOCPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    UOCPIE.final <- sum(UOCPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    UOOPDE.cluster <- c( mean(UOOPPart1[[1]]$DE + UOPart2[[1]]$DE) , mean(UOOPPart1[[2]]$DE + UOPart2[[2]]$DE) )            # Cluster Level DE
    UOOPIE.cluster <- c( mean(UOOPPart1[[1]]$IE + UOPart2[[1]]$IE) , mean(UOOPPart1[[2]]$IE + UOPart2[[2]]$IE) )            # Cluster Level IE
    UOOPDE.final <- sum(UOOPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    UOOPIE.final <- sum(UOOPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    UOUPDE.cluster <- c( mean(UOUPPart1[[1]]$DE + UOPart2[[1]]$DE) , mean(UOUPPart1[[2]]$DE + UOPart2[[2]]$DE) )            # Cluster Level DE
    UOUPIE.cluster <- c( mean(UOUPPart1[[1]]$IE + UOPart2[[1]]$IE) , mean(UOUPPart1[[2]]$IE + UOPart2[[2]]$IE) )            # Cluster Level IE
    UOUPDE.final <- sum(UOUPDE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level DE
    UOUPIE.final <- sum(UOUPIE.cluster*c(Nk[1],Nk[2])/N)                                                            # Population Level IE
    
    
    
    
    OvDE.cluster <- c( mean(OvPart1[[1]]$DE + OvPart2[[1]]$DE) , mean(OvPart1[[2]]$DE + OvPart2[[2]]$DE),
                       mean(OvPart1[[3]]$DE + OvPart2[[3]]$DE) , mean(OvPart1[[4]]$DE + OvPart2[[4]]$DE) )
    OvIE.cluster <- c( mean(OvPart1[[1]]$IE + OvPart2[[1]]$IE) , mean(OvPart1[[2]]$IE + OvPart2[[2]]$IE),
                       mean(OvPart1[[3]]$IE + OvPart2[[3]]$IE) , mean(OvPart1[[4]]$IE + OvPart2[[4]]$IE) )
    OvDE.final <- sum( OvDE.cluster*(wNk/N) )
    OvIE.final <- sum( OvIE.cluster*(wNk/N) )
    
    UdDE.cluster <- mean( c(UdPart1[[1]]$DE + UdPart2[[1]]$DE , UdPart1[[2]]$DE + UdPart2[[2]]$DE) )
    UdIE.cluster <- mean( c(UdPart1[[1]]$IE + UdPart2[[1]]$IE , UdPart1[[2]]$IE + UdPart2[[2]]$IE) )
    UdDE.final <- sum( UdDE.cluster )
    UdIE.final <- sum( UdIE.cluster )
    
    DE.Empvar <- c( var(c(COCPPart1[[1]]$DE + COPart2[[1]]$DE , COCPPart1[[2]]$DE + COPart2[[2]]$DE)),                  # SE estimate of DE
                    var(c(COOPPart1[[1]]$DE + COPart2[[1]]$DE , COOPPart1[[2]]$DE + COPart2[[2]]$DE)),
                    var(c(COUPPart1[[1]]$DE + COPart2[[1]]$DE , COUPPart1[[2]]$DE + COPart2[[2]]$DE)),
                    
                    var(c(OOCPPart1[[1]]$DE + OOPart2[[1]]$DE , OOCPPart1[[2]]$DE + OOPart2[[2]]$DE)), 
                    var(c(OOOPPart1[[1]]$DE + OOPart2[[1]]$DE , OOOPPart1[[2]]$DE + OOPart2[[2]]$DE)),
                    var(c(OOUPPart1[[1]]$DE + OOPart2[[1]]$DE , OOUPPart1[[2]]$DE + OOPart2[[2]]$DE)), 
                    
                    var(c(UOCPPart1[[1]]$DE + UOPart2[[1]]$DE , UOCPPart1[[2]]$DE + UOPart2[[2]]$DE)), 
                    var(c(UOOPPart1[[1]]$DE + UOPart2[[1]]$DE , UOOPPart1[[2]]$DE + UOPart2[[2]]$DE)),
                    var(c(UOUPPart1[[1]]$DE + UOPart2[[1]]$DE , UOUPPart1[[2]]$DE + UOPart2[[2]]$DE)), 
                    
                    var(c(OvPart1[[1]]$DE + OvPart2[[1]]$DE , OvPart1[[2]]$DE + OvPart2[[2]]$DE,
                          OvPart1[[3]]$DE + OvPart2[[3]]$DE , OvPart1[[4]]$DE + OvPart2[[4]]$DE)),
                    
                    var(c(UdPart1[[1]]$DE + UdPart2[[1]]$DE , UdPart1[[2]]$DE + UdPart2[[2]]$DE)) )
    
    IE.Empvar <- c( var(c(COCPPart1[[1]]$IE + COPart2[[1]]$IE , COCPPart1[[2]]$IE + COPart2[[2]]$IE)),                  # SE estimate of IE
                    var(c(COOPPart1[[1]]$IE + COPart2[[1]]$IE , COOPPart1[[2]]$IE + COPart2[[2]]$IE)),
                    var(c(COUPPart1[[1]]$IE + COPart2[[1]]$IE , COUPPart1[[2]]$IE + COPart2[[2]]$IE)),
                    
                    var(c(OOCPPart1[[1]]$IE + OOPart2[[1]]$IE , OOCPPart1[[2]]$IE + OOPart2[[2]]$IE)), 
                    var(c(OOOPPart1[[1]]$IE + OOPart2[[1]]$IE , OOOPPart1[[2]]$IE + OOPart2[[2]]$IE)),
                    var(c(OOUPPart1[[1]]$IE + OOPart2[[1]]$IE , OOUPPart1[[2]]$IE + OOPart2[[2]]$IE)), 
                    
                    var(c(UOCPPart1[[1]]$IE + UOPart2[[1]]$IE , UOCPPart1[[2]]$IE + UOPart2[[2]]$IE)), 
                    var(c(UOOPPart1[[1]]$IE + UOPart2[[1]]$IE , UOOPPart1[[2]]$IE + UOPart2[[2]]$IE)),
                    var(c(UOUPPart1[[1]]$IE + UOPart2[[1]]$IE , UOUPPart1[[2]]$IE + UOPart2[[2]]$IE)), 
                    
                    var(c(OvPart1[[1]]$IE + OvPart2[[1]]$IE , OvPart1[[2]]$IE + OvPart2[[2]]$IE,
                          OvPart1[[3]]$IE + OvPart2[[3]]$IE , OvPart1[[4]]$IE + OvPart2[[4]]$IE)),
                    
                    var(c(UdPart1[[1]]$IE + UdPart2[[1]]$IE , UdPart1[[2]]$IE + UdPart2[[2]]$IE)) )
    
    Result[iteration,1:44] <- c( COCPDE.final, COOPDE.final, COUPDE.final,
                                 OOCPDE.final, OOOPDE.final, OOUPDE.final,
                                 UOCPDE.final, UOOPDE.final, UOUPDE.final,
                                 OvDE.final, UdDE.final,
                                 COCPIE.final, COOPIE.final, COUPIE.final,
                                 OOCPIE.final, OOOPIE.final, OOUPIE.final,
                                 UOCPIE.final, UOOPIE.final, UOUPIE.final,
                                 OvIE.final, UdIE.final,
                                 DE.Empvar, IE.Empvar)
    
    Result[iteration,45:66] <- 
        as.numeric( Result[iteration,1:22] - qnorm(0.975)*sqrt(Result[iteration,23:44])/sqrt(N) - rep(c(TrueDE, TrueIE),each=11) <= 0) *  
        as.numeric( Result[iteration,1:22] + qnorm(0.975)*sqrt(Result[iteration,23:44])/sqrt(N) - rep(c(TrueDE, TrueIE),each=11) >= 0) 
    
    Result[iteration,1:22] <- Result[iteration,1:22] - rep(c(TrueDE, TrueIE),each=11)
    
    print(iteration)
    
}

Result <- data.frame(Result)
colnames(Result) <- c("B_True_DE","B_MisOR_DE","B_MisPS_DE","B_MisBoth_DE","B_Over_DE","B_Under_DE",
                      "B_True_IE","B_MisOR_IE","B_MisPS_IE","B_MisBoth_IE","B_Over_IE","B_Under_IE",
                      "S_True_DE","S_MisOR_DE","S_MisPS_DE","S_MisBoth_DE","S_Over_DE","S_Under_DE",
                      "S_True_IE","S_MisOR_IE","S_MisPS_IE","S_MisBoth_IE","S_Over_IE","S_Under_IE",
                      "C_True_DE","C_MisOR_DE","C_MisPS_DE","C_MisBoth_DE","C_Over_DE","C_Under_DE",
                      "C_True_IE","C_MisOR_IE","C_MisPS_IE","C_MisBoth_IE","C_Over_IE","C_Under_IE")

write.csv(Result,"RESULT16.csv",row.names=FALSE)                                   # Save Result



############################
# Report the Result
############################

Result <- read.csv("RESULT16.csv")                                         # Load the Provided Datasets

position <- c(1,3,7,9,11,10)

N <- 2000                                                                        # Number of Clusters

cbind( apply(Result,2,mean)[position],
       apply(Result,2,sd)[position]/sqrt(N),
       apply(Result,2,mean)[44+position],
       apply(Result,2,mean)[11+position],
       apply(Result,2,sd)[11+position]/sqrt(N),
       apply(Result,2,mean)[55+position] )


