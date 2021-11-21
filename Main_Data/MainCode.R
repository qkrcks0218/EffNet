############################################
# Last update : October 02 2021
# Code for the data analysis in Section 5
############################################

##################################################
# Load Packages and Functions
##################################################

library(readstata13)
library(SuperLearner)
library(caret)
library(nnet)
library(glmnet)
library(earth)
library(gam)
library(gbm)
library(xgboost)     
library(kernlab)
library(polspline)
library(ranger)
source("MainFunction.R")

##################################################
# Data Preparation: Download the data from the link below
# https://www.openicpsr.org/openicpsr/project/113783/version/V1/view
##################################################

RawData <- read.dta13("Public_Data_AEJApp_2010-0132.dta")

################################
# Calculate Propensity Score
################################

SC <- as.numeric(RawData$suba==0)                       # Student in San Cristobal. These students were eligible for the lottery.
SC.eligible <- SC                                       # Eligible student in San Cristobal.
SUBA <- RawData$suba                                    # Student in Suba. Students under grade 9 were not eligible for the lottery.
SUBA.eligible <- SUBA * as.numeric(RawData$grade>=9)    # Eligible student in Suba. 
Total.SC <- sum(SC.eligible)                            # 10,907. Agrees with  Barrera-Osorio et al. (2011)
Total.SUBA <- sum(SUBA.eligible)                        # 2,526. Agrees with  Barrera-Osorio et al. (2011)

Trt.SC <- sum(SC.eligible * RawData$treatment)          # 6,851. Total number of students received treatment in San Cristobal.
Trt.SUBA <- sum(SUBA.eligible * RawData$treatment)      # 1,133. Total number of students received treatment in San Cristobal.

# Strata used in the randomization

SC.Pub.M.06 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==06)
SC.Pub.M.07 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==07)
SC.Pub.M.08 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==08)
SC.Pub.M.09 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==09)
SC.Pub.M.10 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==10)
SC.Pub.M.11 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==11)

SC.Pub.F.06 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==06)
SC.Pub.F.07 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==07)
SC.Pub.F.08 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==08)
SC.Pub.F.09 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==09)
SC.Pub.F.10 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==10)
SC.Pub.F.11 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==11)

SC.Pri.M.06 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==06)
SC.Pri.M.07 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==07)
SC.Pri.M.08 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==08)
SC.Pri.M.09 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==09)
SC.Pri.M.10 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==10)
SC.Pri.M.11 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==11)

SC.Pri.F.06 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==06)
SC.Pri.F.07 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==07)
SC.Pri.F.08 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==08)
SC.Pri.F.09 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==09)
SC.Pri.F.10 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==10)
SC.Pri.F.11 <- as.numeric(RawData$suba==0 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==11)

Suba.Pub.M.09 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==09)
Suba.Pub.M.10 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==10)
Suba.Pub.M.11 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code >1 & RawData$grade==11)

Suba.Pub.F.09 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==09)
Suba.Pub.F.10 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==10)
Suba.Pub.F.11 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code >1 & RawData$grade==11)

Suba.Pri.M.09 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==09)
Suba.Pri.M.10 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==10)
Suba.Pri.M.11 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="M" & RawData$school_code==1 & RawData$grade==11)

Suba.Pri.F.09 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==09)
Suba.Pri.F.10 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==10)
Suba.Pri.F.11 <- as.numeric(RawData$suba==1 & RawData$r_be_gene=="F" & RawData$school_code==1 & RawData$grade==11)

Randomization <- cbind(SC.Pub.M.06,SC.Pub.M.07,SC.Pub.M.08,SC.Pub.M.09,SC.Pub.M.10,SC.Pub.M.11,
                       SC.Pub.F.06,SC.Pub.F.07,SC.Pub.F.08,SC.Pub.F.09,SC.Pub.F.10,SC.Pub.F.11,
                       SC.Pri.M.06,SC.Pri.M.07,SC.Pri.M.08,SC.Pri.M.09,SC.Pri.M.10,SC.Pri.M.11,
                       SC.Pri.F.06,SC.Pri.F.07,SC.Pri.F.08,SC.Pri.F.09,SC.Pri.F.10,SC.Pri.F.11,
                       Suba.Pub.M.09,Suba.Pub.M.10,Suba.Pub.M.11,
                       Suba.Pub.F.09,Suba.Pub.F.10,Suba.Pub.F.11,
                       Suba.Pri.M.09,Suba.Pri.M.10,Suba.Pri.M.11,
                       Suba.Pri.F.09,Suba.Pri.F.10,Suba.Pri.F.11)

RP1 <- apply(Randomization*RawData$treatment,2,sum)
RP2 <- apply(Randomization,2,sum)
RP1 <- apply(Randomization*matrix(RP1,dim(Randomization)[1],dim(Randomization)[2],byrow=T),1,sum)
RP2 <- apply(Randomization*matrix(RP2,dim(Randomization)[1],dim(Randomization)[2],byrow=T),1,sum)
RawData <- cbind(RawData,Randomization,RP1,RP2)

################################
# Data Cleaning
################################

Data <- RawData[RawData$survey_selected==1&             # Students who enrolled in one of the 68 schools selected for surveying and verified attendance.
                    (SC.eligible+SUBA.eligible)==1,]    # and who were eligible for the lottery. 7,569 students.

################################
# Table
################################

Data.SUBA <- Data[Data$suba==1,]
Data.SC <- Data[Data$suba==0,]

################################
# Cluster ID (Household) Assignment
# Some cluster ID is missing and
# some students under the same cluster ID have different cluster-level variables.
# We remove clusters of which students report different cluster-level variables
################################

Cluster.cov.list <- c("s_edadhead",                     # Age of household head
                      "s_single",                       # Single parent household
                      "s_tpersona",                     # Number of people in household
                      "s_puntaje",                      # SISBEN score
                      "s_ingtotal",                     # Household income
                      "suba")                           # Location=Suba

Data <- Data[!is.na(Data$fu_nim_hogar),]                # Choose students with known cluster ID
Data <- Data[order(Data$fu_nim_hogar),]                 # Sort based on cluster ID

CID <- unique(Data$fu_nim_hogar)                        # Cluster ID

CX <- Data[,Cluster.cov.list]
valid <- rep(1,dim(Data)[1])

for(ii in 1:length(CID)){                               # valid = 0 if cluster-level variables are inconsistent across siblings
    INDEX <- which(Data$fu_nim_hogar==CID[ii])
    if( length(INDEX)>1&sum(apply(CX[INDEX,],2,var))>1 ){
        valid[INDEX] <- 0
    }
}

Data <- Data[valid==1,]                                 # Choose valid data

Data$fu_nim_hogar <-                                    # Re-assign cluster ID
    rep(1:length(unique(Data$fu_nim_hogar)),table(Data$fu_nim_hogar))

Data$school_code <- as.numeric(Data$school_code==1)      # Private school

################################
# Covariate Selection
################################

RDN <- colnames(Data)[51:88]
Covariate.list <- c("s_age_sorteo",                     # Age
                    "grade",                            # Grade
                    "s_sexo",                           # Gender
                    "s_edadhead",                       # Age of household head
                    "s_single",                         # Single parent household
                    "s_tpersona",                       # Number of people in household
                    "s_puntaje",                        # SISBEN score
                    "s_ingtotal",                       # Household income
                    "suba",                             # Location=Suba
                    "fu_nim_hogar",                     # Household ID
                    RDN)                                # Randomization

# Data <- Data[,c("m_enrolled","treatment",               # Enrollment, Treatment
#                 Covariate.list)]
# Data <- Data[,c("at_msamean","treatment",               # Monitored attendance, Treatment
#                 Covariate.list)]
Data <- Data[,c("fu_self_attendance","treatment",               # Self attendance, Treatment
                Covariate.list)]
Data <- Data[apply(is.na(Data),1,sum)==0,]              # Choose Observations that have full information.
Data <- Data[order(Data$fu_nim_hogar),]                 # Sort based on household ID

colnames(Data) <- c("Y","A",                            # Rename
                    paste("Ind.",c("Age","Grade","Gender"),sep=""),
                    paste("HH.",c("Age","FamilyType","NumPpl","SISBEN","Income","Loc","ID"),sep=""),
                    RDN)


M <- as.numeric( table( Data$HH.ID ) )
Data$M <- rep(M,M)                                      # Cluster size

write.csv(Data[,-c(13:50)],"Clean_Data_SelfAttendance.csv",row.names=FALSE)
write.csv(Data[,c(13:50)],"Clean_Data_SelfAttendance_Randomization.csv",row.names=FALSE)


################################
# Load Data
################################

Data <- read.csv("Clean_Data_SelfAttendance.csv")

################################
# Superlearner Hyperparameter Grids
################################

SL.hpara <- list()
SL.hpara$SLL <- 1 #:9
# Superlearner basic learning algorithms: 
# 1: GLM
# 2: lasso/ridge
# 3: earth
# 4: GAM
# 5: xgboost
# 6: polynomial spline
# 7: random forest
# 8: gbm
# 9: 1-layer MLP
SL.hpara$MLPL <- c(2,4)                 # Hyperparameter for 9: 1-layer MLP
SL.hpara$MLPdecay <- 10^c(-1,-3,-5)     # Hyperparameter for 9: 1-layer MLP
SL.hpara$MTRY <- c(3,6)                 # Hyperparameter for 7: random forest
ORtype <- gaussian()                    # Specify continuous outcome

################################
# Type Definition: Product of M * Location
################################

TypeD <- list()
TypeD[[1]] <- Data[Data$M==1&Data$HH.Loc==0,]                   # Size 1 & San Cristobal
TypeD[[2]] <- Data[Data$M==1&Data$HH.Loc==1,]                   # Size 1 & Suba
TypeD[[3]] <- Data[Data$M==2&Data$HH.Loc==0,]                   # Size 2 & San Cristobal
TypeD[[4]] <- Data[Data$M==2&Data$HH.Loc==1,]                   # Size 2 & Suba
TypeD[[5]] <- Data[Data$M>=3&Data$HH.Loc==0,]                   # Size>=3& San Cristobal
TypeD[[6]] <- Data[Data$M>=3&Data$HH.Loc==1,]                   # Size>=3& Suba


N <- M <- Index <- C <- list()
for(tt in 1:6){
    C[[tt]] <- unique(TypeD[[tt]]$HH.ID)                              # Cluster list in tt-th type
    N[[tt]] <- length(C[[tt]])                                        # Number of clusters
    Index[[tt]] <- matrix(0, length(C[[tt]]), 2)                      # Index of each cluster
    for(ii in 1:dim(Index[[tt]])[1]){
        Index[[tt]][ii,1] <- min( which(TypeD[[tt]]$HH.ID==C[[tt]][ii]) )
        Index[[tt]][ii,2] <- max( which(TypeD[[tt]]$HH.ID==C[[tt]][ii]) )
    }
    M[[tt]] <- Index[[tt]][,2]-Index[[tt]][,1]+1                      # Cluster size
}



for(jj in 3:6){                                                 # Aggregate others' A and X
    Others <- aggregate(.~HH.ID,data=TypeD[[jj]],FUN="sum")[,c(3:6,1,13)]
    Others.Expand <- matrix(0,dim(TypeD[[jj]])[1],4)
    for(ii in 1:dim(Others)[1]){
        Others.Expand[Index[[jj]][ii,1]:Index[[jj]][ii,2],] <- matrix(as.numeric(Others[ii,1:4]),M[[jj]][ii],4,byrow=T)
    }
    Others.Summary <- as.matrix((Others.Expand - as.matrix(TypeD[[jj]][,c(2:5)]))/matrix(TypeD[[jj]]$M-1,dim(TypeD[[jj]])[1],4))
    Others.Summary <- data.frame(Others.Summary)
    colnames(Others.Summary) <- paste("Other.",c("A","Age","Grade","Gender"),sep="")
    # Others.Summary$Other.A.Ego.A <- Others.Summary$Other.A*TypeD[[jj]]$A
    TypeD[[jj]] <- cbind(TypeD[[jj]][,c(1:5)],Others.Summary,TypeD[[jj]][,c(6:13)])
}

for(jj in 1:2){ TypeD[[jj]] <- TypeD[[jj]][,-c(which(colnames(TypeD[[jj]])=="HH.Loc"),
                                               which(colnames(TypeD[[jj]])=="M"))] }            # Remove Location, M
for(jj in 3:4){ TypeD[[jj]] <- TypeD[[jj]][,-c(which(colnames(TypeD[[jj]])=="HH.Loc"),
                                               which(colnames(TypeD[[jj]])=="M"))] }            # Remove Location, M
for(jj in 5:6){ TypeD[[jj]] <- TypeD[[jj]][,-c(which(colnames(TypeD[[jj]])=="HH.Loc"))] }       # Remove Location

################################
# Propensity Score
################################

Randomization <- read.csv("Clean_Data_SelfAttendance_Randomization.csv")
TypeR <- list()
TypeR[[1]] <- Randomization[Data$M==1&Data$HH.Loc==0,]                   # Size 1 & San Cristobal
TypeR[[2]] <- Randomization[Data$M==1&Data$HH.Loc==1,]                   # Size 1 & Suba
TypeR[[3]] <- Randomization[Data$M==2&Data$HH.Loc==0,]                   # Size 2 & San Cristobal
TypeR[[4]] <- Randomization[Data$M==2&Data$HH.Loc==1,]                   # Size 2 & Suba
TypeR[[5]] <- Randomization[Data$M>=3&Data$HH.Loc==0,]                   # Size>=3& San Cristobal
TypeR[[6]] <- Randomization[Data$M>=3&Data$HH.Loc==1,]                   # Size>=3& Suba

PS <- list()
PS[[1]] <- (TypeR[[1]][,37]/TypeR[[1]][,38])^TypeD[[1]]$A*(1-TypeR[[1]][,37]/TypeR[[1]][,38])^(1-TypeD[[1]]$A)
PS[[2]] <- (TypeR[[2]][,37]/TypeR[[2]][,38])^TypeD[[2]]$A*(1-TypeR[[2]][,37]/TypeR[[2]][,38])^(1-TypeD[[2]]$A)

for(tt in 3:4){
    PS[[tt]] <- rep(0,dim(TypeD[[tt]])[1])                         # Propensity score of observing treatment assignment vectors 
    for(ii in 1:dim(Index[[tt]])[1]){
        A.temp <- TypeD[[tt]]$A[Index[[tt]][ii,1]:Index[[tt]][ii,2]]
        R.temp <- TypeR[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2],]
        if( sum(apply(R.temp,2,sum)[-c(37,38)]>1)==0 ){
            PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- prod( (R.temp[,37]/R.temp[,38])^A.temp*(1-R.temp[,37]/R.temp[,38])^(1-A.temp) )
        } else {
            if(sum(A.temp)==2){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- (R.temp[1,37]*(R.temp[2,37]-1))/R.temp[1,38]/(R.temp[1,38]-1)
            } else if (sum(A.temp)==1){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- (R.temp[1,37]*(R.temp[1,38]-R.temp[2,37]))/R.temp[1,38]/(R.temp[1,38]-1)
            } else if (sum(A.temp)==0){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- (R.temp[1,38]-R.temp[1,37])*(R.temp[1,38]-R.temp[2,37]-1)/R.temp[1,38]/(R.temp[1,38]-1)
            }
        }
        
    }
}


for(tt in 5){
    PS[[tt]] <- rep(0,dim(TypeD[[tt]])[1])                         # Propensity score of observing treatment assignment vectors 
    for(ii in 1:dim(Index[[tt]])[1]){
        A.temp <- TypeD[[tt]]$A[Index[[tt]][ii,1]:Index[[tt]][ii,2]]
        R.temp <- TypeR[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2],]
        if( sum(apply(R.temp,2,sum)[-c(37,38)]>1)==0 ){
            PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- prod( (R.temp[,37]/R.temp[,38])^A.temp*(1-R.temp[,37]/R.temp[,38])^(1-A.temp) )
        } else {
            pos.2 <- which(R.temp[,which(apply(R.temp,2,sum)[-c(37,38)]==2)]==1)
            pos.1 <- (1:length(A.temp))[-pos.2]
            p1 <- prod( (R.temp[pos.1,37]/R.temp[pos.1,38])^A.temp[pos.1]*(1-R.temp[pos.1,37]/R.temp[pos.1,38])^(1-A.temp[pos.1]) )
            A.temp <- A.temp[pos.2]
            R.temp <- R.temp[pos.2,]
            if(sum(A.temp)==2){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- p1*(R.temp[1,37]*(R.temp[2,37]-1))/R.temp[1,38]/(R.temp[1,38]-1)
            } else if (sum(A.temp)==1){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- p1*(R.temp[1,37]*(R.temp[1,38]-R.temp[2,37]))/R.temp[1,38]/(R.temp[1,38]-1)
            } else if (sum(A.temp)==0){
                PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- p1*(R.temp[1,38]-R.temp[1,37])*(R.temp[1,38]-R.temp[2,37]-1)/R.temp[1,38]/(R.temp[1,38]-1)
            }
        }
    }
    
}


for(tt in 6){
    PS[[tt]] <- rep(0,dim(TypeD[[tt]])[1])                         # Propensity score of observing treatment assignment vectors 
    for(ii in 1:dim(Index[[tt]])[1]){
        A.temp <- TypeD[[tt]]$A[Index[[tt]][ii,1]:Index[[tt]][ii,2]]
        R.temp <- TypeR[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2],]
        if( sum(apply(R.temp,2,sum)[-c(37,38)]>1)==0 ){
            PS[[tt]][Index[[tt]][ii,1]:Index[[tt]][ii,2]] <- prod( (R.temp[,37]/R.temp[,38])^A.temp*(1-R.temp[,37]/R.temp[,38])^(1-A.temp) )
        } else {
            print(c(tt,ii))
        }
    }
}

#################################################
# ML Analysis : Only use type 3,4,5
#################################################

for(SS.ITERATION in 1:100){         # Repeat sampls splitting 100 times
    
    
    SS <- list()                    # Cluster-level split sample indices
    SS[[3]] <- SS[[4]] <- SS[[5]] <- list()
    SSI <- list()                   # Individual-level split sample indices
    SSI[[3]] <- SSI[[4]] <- SSI[[5]] <- list()
    
    for(ii in 3:5){
        SS[[ii]]$MS <- sort(sample(N[[ii]],round(N[[ii]]/2)))       # Sample 1
        SS[[ii]]$AS <- (1:N[[ii]])[-SS[[ii]]$MS]                    # Sample 2
        if(ii<=4) {
            SSI[[ii]]$MS <- NULL
            SSI[[ii]]$AS <- NULL
            for(jj in 1:N[[ii]]){
                if( sum(jj==SS[[ii]]$MS)>0 ){
                    SSI[[ii]]$MS <- c(SSI[[ii]]$MS,Index[[ii]][jj,1]:Index[[ii]][jj,2])
                } else {
                    SSI[[ii]]$AS <- c(SSI[[ii]]$AS,Index[[ii]][jj,1]:Index[[ii]][jj,2])
                }
            }
        } else {
            S2 <- which(Index[[5]][,2]-Index[[5]][,1]==2)
            S3 <- which(Index[[5]][,2]-Index[[5]][,1]>2)
            SS[[ii]]$MS <- sort(c(sample(S2,length(S2)/2), sample(S3,round(length(S3)/2))))
            SS[[ii]]$AS <- (1:N[[ii]])[-SS[[ii]]$MS]
            SSI[[ii]]$MS <- NULL
            SSI[[ii]]$AS <- NULL
            for(jj in (1:N[[ii]])){
                if( sum(jj==SS[[ii]]$MS)>0 ){
                    SSI[[ii]]$MS <- c(SSI[[ii]]$MS,Index[[ii]][jj,1]:Index[[ii]][jj,2])
                } else {
                    SSI[[ii]]$AS <- c(SSI[[ii]]$AS,Index[[ii]][jj,1]:Index[[ii]][jj,2])
                }
            }
        }
    }
    
    
    OR.fit <- OR.est <- DB <- list()
    OR.fit[[3]] <- OR.fit[[4]] <- OR.fit[[5]] <- list()             # Estimated OR function
    OR.est[[3]] <- OR.est[[4]] <- OR.est[[5]] <- list()             # Evaluated OR estimate
    DB <- RESULT <- list()                                          # Doubly robust term and result
    
    for(ii in 3:5){
        
        TempData <- TypeD[[ii]]
        if(ii<=4){
            TempData <- TempData[,-c(which(colnames(TempData)=="HH.Loc"),which(colnames(TempData)=="HH.ID"),which(colnames(TempData)=="M"))]
        } else {
            TempData <- TempData[,-c(which(colnames(TempData)=="HH.Loc"),which(colnames(TempData)=="HH.ID"))]
            TempData$Other.A.C1 <- as.numeric(TempData$Other.A==0)
            TempData$Other.A.C2 <- as.numeric(TempData$Other.A>0 & TempData$Other.A<=0.5)
            TempData$Other.A.C3 <- as.numeric(TempData$Other.A>0.5 & TempData$Other.A<=0.75)
            TempData <- TempData[,-c(which(colnames(TempData)=="Other.A"))]
        }
        
        AmatforType <- makingAmat(max(M[[ii]]))
        
        TempData.ReplaceA <- list()
        for(att in 1:dim(AmatforType)[1] ){
            TempData.ReplaceA[[att]] <- TempData
        }
        for(jj in 1:N[[ii]]){
            pos <- Index[[ii]][jj,1]:Index[[ii]][jj,2]
            for(att in 1:(2^M[[ii]][jj]) ){
                TempData.ReplaceA[[att]]$A[pos] <-  makingAmat(M[[ii]][jj])[att,]
                if(3<=ii & ii <=4){
                    for(ll in 1:length(pos)){
                        TempData.ReplaceA[[att]]$Other.A[pos[ll]] <- mean(TempData.ReplaceA[[att]]$A[pos[-ll]])
                    }
                } else if (ii==5){
                    OtherProp <- mean(TempData.ReplaceA[[att]]$A[pos[-ll]])
                    TempData.ReplaceA[[att]]$Other.A.C1[pos[ll]] <- as.numeric(OtherProp==0)
                    TempData.ReplaceA[[att]]$Other.A.C2[pos[ll]] <- as.numeric(OtherProp>0 & OtherProp<=0.5)
                    TempData.ReplaceA[[att]]$Other.A.C3[pos[ll]] <- as.numeric(OtherProp>0.5 & OtherProp<=0.75)
                    
                }
            }
            if(M[[ii]][jj] < max(M[[ii]])){
                for(att in (2^M[[ii]][jj]+1):(dim(AmatforType)[1]) ){
                    TempData.ReplaceA[[att]]$A[pos] <-  NA
                    if(3<= ii & ii<=4){ TempData.ReplaceA[[att]]$Other.A[pos] <-  NA }
                    if(ii==5){ TempData.ReplaceA[[att]]$Other.A.C1[pos] <-  NA
                    TempData.ReplaceA[[att]]$Other.A.C2[pos] <-  NA
                    TempData.ReplaceA[[att]]$Other.A.C3[pos] <-  NA}
                }
            }
        }
        
        
        SL.hpara$NMN <- round(0.05*N[[ii]])  # Hyperparameter for 8: gbm
        pos.X <- 2:dim(TempData)[2]          # Position of X
        
        ############### Outcome Model Estimation via Superlearner
        
        OR.fit[[ii]]$MS <- MySL(TempData[SSI[[ii]]$MS,],locY=1,locX=pos.X,Ydist=ORtype,
                                SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, 
                                NMN=SL.hpara$NMN, MLPdecay=SL.hpara$MLPdecay)
        OR.fit[[ii]]$AS <- MySL(TempData[SSI[[ii]]$AS,],locY=1,locX=pos.X,Ydist=ORtype,
                                SL.list=SL.hpara$SLL, MTRY=SL.hpara$MTRY, MLPL=SL.hpara$MLPL, 
                                NMN=SL.hpara$NMN, MLPdecay=SL.hpara$MLPdecay)
        
        ############### Outcome Model Estimates
        
        OR.est[[ii]] <- rep(0,N[[ii]])
        OR.est[[ii]][SSI[[ii]]$AS] <- predict(OR.fit[[ii]]$MS,TempData[SSI[[ii]]$AS,pos.X])$pred
        OR.est[[ii]][SSI[[ii]]$MS] <- predict(OR.fit[[ii]]$AS,TempData[SSI[[ii]]$MS,pos.X])$pred
        
        OR.est.newA <- list()
        for(att in 1:dim(AmatforType)[1] ){
            valid <- which(!is.na(TempData.ReplaceA[[att]]$A))
            OR.est.newA[[att]] <- rep(NA,sum(M[[ii]]))
            if ( length( intersect(SSI[[ii]]$AS,valid) ) > 0 ){
                OR.est.newA[[att]][intersect(SSI[[ii]]$AS,valid)] <- predict(OR.fit[[ii]]$MS,TempData.ReplaceA[[att]][intersect(SSI[[ii]]$AS,valid),pos.X])$pred
            }
            if ( length( intersect(SSI[[ii]]$MS,valid) ) > 0 ){
                OR.est.newA[[att]][intersect(SSI[[ii]]$MS,valid)] <- predict(OR.fit[[ii]]$AS,TempData.ReplaceA[[att]][intersect(SSI[[ii]]$MS,valid),pos.X])$pred
            }
        }
        
        ############### Weighted Residual Term 
        
        DB[[ii]] <- (TempData$Y-OR.est[[ii]])/PS[[ii]]
        
        ############### Full composition of estimated residual/outcome regression
        
        RESULT[[ii]] <- cbind(DB[[ii]],OR.est[[ii]])
        
        for(att in 1:dim(AmatforType)[1] ){
            RESULT[[ii]] <- cbind(RESULT[[ii]],OR.est.newA[[att]])   
        }
        
        if(ii<=2){
            colnames(RESULT[[ii]]) <- c("DB","OR","OR0","OR1")
        } else if (ii<=4) {
            colnames(RESULT[[ii]]) <- c("DB","OR","OR00","OR01","OR10","OR11")
        } else if (ii==5) { 
            colnames(RESULT[[ii]]) <- c("DB","OR",paste("OR",paste(AmatforType[,1],AmatforType[,2],AmatforType[,3],AmatforType[,4],AmatforType[,5],sep=""),sep=""))
        }
        
    }
    
    # Save Data: this will make csv files in Estimate Folder

    # for(ii in 3:5){
    #     write.csv(RESULT[[ii]],sprintf("Estimate/RESULT_Disc_B%0.4d_T%0.1d.csv",SS.ITERATION,ii),row.names=FALSE)
    # }
    
}

#################################################
# Summarizing Results
#################################################

RESULT <- list()
for(ii in 1:100){
    RESULT[[ii]] <- list()
    RESULT[[ii]][[3]] <- read.csv(sprintf("Estimate/RESULT_Disc_B%0.4d_T3.csv",ii))
    RESULT[[ii]][[4]] <- read.csv(sprintf("Estimate/RESULT_Disc_B%0.4d_T4.csv",ii))
    RESULT[[ii]][[5]] <- read.csv(sprintf("Estimate/RESULT_Disc_B%0.4d_T5.csv",ii))
}

#################################################
# Obtain Effect, SE, Wald statistic over [0,1]^2
#################################################

Eff.IPW.DE <- Eff.IPW.IE <- Eff.ML.DE <- Eff.ML.IE <- matrix(0,21,21)
SE.IPW.DE <- SE.IPW.IE <- SE.ML.DE <- SE.ML.IE <- matrix(0,21,21)
Wald.IPW.DE <- Wald.IPW.IE <- Wald.ML.DE <- Wald.ML.IE <- matrix(0,21,21)

g1 <- (0:20)/20     # grid for alpha_{SC}
g2 <- (0:20)/20     # grid for alpha_{Suba}

for(index1 in 1:length(g1)){
    
    for(index2 in 1:length(g2)){
        
        #################################################
        # IPW Analysis
        #################################################
        
        a1 <- g1[index1]
        a2 <- g2[index2]
    
        # Counterfactual parameters
        alpha.grid <- c(NA,NA,a1,a2,a1,NA)
        alpha1.grid <- c(NA,NA,a1,a2,a1,NA)
        alpha2.grid <- c(NA,NA,Trt.SC/Total.SC,Trt.SUBA/Total.SUBA,Trt.SC/Total.SC,NA)
        
        IPW.IF.DE <- list()
        IPW.IF.IE <- list()
        
        Wvec.DE <- list()
        Wvec.DE[[3]] <- Wvec.DE[[4]] <- Wvec.DE[[5]] <- list()
        Wvec.IE <- list()
        Wvec.IE[[3]] <- Wvec.IE[[4]] <- Wvec.IE[[5]] <- list()
        
        
        # IPW Estimator
        for(ii in 3:5){
            TempData <- TypeD[[ii]]
            
            Wvec.DE[[ii]] <- rep(0,N[[ii]])
            Wvec.IE[[ii]] <- rep(0,N[[ii]])
            
            for(i in 1:N[[ii]]){
                Mk <- Index[[ii]][i,2]-Index[[ii]][i,1]+1
                WMat <- DEIEmat(Mk,alpha=alpha.grid[ii], alpha1=alpha1.grid[ii], alpha2=alpha2.grid[ii])
                Avec <- TempData$A[Index[[ii]][i,1]:Index[[ii]][i,2]]
                Wvec.DE[[ii]][Index[[ii]][i,1]:Index[[ii]][i,2]] <- WMat$DE[sum(Avec*2^((Mk-1):0))+1,]
                Wvec.IE[[ii]][Index[[ii]][i,1]:Index[[ii]][i,2]] <- WMat$IE[sum(Avec*2^((Mk-1):0))+1,]
            }
            
            IPW.IF.DE[[ii]] <- aggregate( Wvec.DE[[ii]]*TempData$Y/PS[[ii]]~TempData$HH.ID, FUN="sum")[,2]
            IPW.IF.IE[[ii]] <- aggregate( Wvec.IE[[ii]]*TempData$Y/PS[[ii]]~TempData$HH.ID, FUN="sum")[,2]
            
            
        }
        
        IPW.IF.DE <- c(IPW.IF.DE[[3]],IPW.IF.DE[[4]],IPW.IF.DE[[5]])
        IPW.IF.IE <- c(IPW.IF.IE[[3]],IPW.IF.IE[[4]],IPW.IF.IE[[5]])
        
        IPW.Est.DE <- c( mean(IPW.IF.DE[[3]]), mean(IPW.IF.DE[[4]]), 
                             mean(IPW.IF.DE[[5]]), mean(IPW.IF.DE) )
        IPW.Est.IE <- c( mean(IPW.IF.IE[[3]]), mean(IPW.IF.IE[[4]]), 
                             mean(IPW.IF.IE[[5]]), mean(IPW.IF.IE) )
        IPW.SE.DE <- c( sd(IPW.IF.DE[[3]])/length(IPW.IF.DE[[3]]),
                            sd(IPW.IF.DE[[4]])/length(IPW.IF.DE[[4]]),
                            sd(IPW.IF.DE[[5]])/length(IPW.IF.DE[[5]]),
                            sd(IPW.IF.DE)/sqrt(length(IPW.IF.DE)) )
        IPW.SE.IE <- c( sd(IPW.IF.IE[[3]])/length(IPW.IF.IE[[3]]),
                            sd(IPW.IF.IE[[4]])/length(IPW.IF.IE[[4]]),
                            sd(IPW.IF.IE[[5]])/length(IPW.IF.IE[[5]]),
                            sd(IPW.IF.IE)/sqrt(length(IPW.IF.IE)) )
        
        ################################
        # ML Method
        ################################
        
        ML.Est.DE.Mat <- ML.Est.IE.Mat <- ML.SE.DE.Mat <- ML.SE.IE.Mat <- matrix(0,100,4)
        
        for(ii in 1:100){
            ML.IF.DE <- ML.IF.IE <- list()
            ML.EST.DE <- ML.EST.IE <- ML.SE.DE <- ML.SE.IE <- rep(0,3)
            
            # Type 3,4, and 5: Size >= 2
            for(tt in 3:5){
                ML.IF.DE[[tt]] <- ML.IF.IE[[tt]] <- rep(0,N[[tt]])
                for(jj in 1:N[[tt]]){
                    pos <- Index[[tt]][jj,1]:Index[[tt]][jj,2]
                    DImat <- DEIEmat(M[[tt]][jj],alpha=alpha.grid[tt], alpha1=alpha1.grid[tt], alpha2=alpha2.grid[tt])
                    Apos <- sum( TypeD[[tt]]$A[pos]*2^(seq(M[[tt]][jj]-1,0,by=-1)) ) + 1
                    ML.IF.DE[[tt]][jj] <- sum(DImat$DE[Apos,]*RESULT[[ii]][[tt]][pos,1]) + sum( t( RESULT[[ii]][[tt]][pos,2+1:(2^M[[tt]][jj])] )*DImat$DE )
                    ML.IF.IE[[tt]][jj] <- sum(DImat$IE[Apos,]*RESULT[[ii]][[tt]][pos,1]) + sum( t( RESULT[[ii]][[tt]][pos,2+1:(2^M[[tt]][jj])] )*DImat$IE )
                }
            }
            
            ML.IF.DE <- c(ML.IF.DE[[3]],ML.IF.DE[[4]],ML.IF.DE[[5]])
            ML.IF.IE <- c(ML.IF.IE[[3]],ML.IF.IE[[4]],ML.IF.IE[[5]])
            
            ML.Est.DE.Mat[ii,] <- c( mean(ML.IF.DE[[3]]), mean(ML.IF.DE[[4]]), 
                                     mean(ML.IF.DE[[5]]), mean(ML.IF.DE) )
            ML.Est.IE.Mat[ii,] <- c( mean(ML.IF.IE[[3]]), mean(ML.IF.IE[[4]]), 
                                     mean(ML.IF.IE[[5]]), mean(ML.IF.IE) )
            ML.SE.DE.Mat[ii,] <- c( sd(ML.IF.DE[[3]])/length(ML.IF.DE[[3]]),
                                    sd(ML.IF.DE[[4]])/length(ML.IF.DE[[4]]),
                                    sd(ML.IF.DE[[5]])/length(ML.IF.DE[[5]]),
                                    sd(ML.IF.DE)/sqrt(length(ML.IF.DE)) )
            ML.SE.IE.Mat[ii,] <- c( sd(ML.IF.IE[[3]])/length(ML.IF.IE[[3]]),
                                    sd(ML.IF.IE[[4]])/length(ML.IF.IE[[4]]),
                                    sd(ML.IF.IE[[5]])/length(ML.IF.IE[[5]]),
                                    sd(ML.IF.IE)/sqrt(length(ML.IF.IE)) )
            print(c(index1,index2,ii))
        }
        
        ML.Est.DE <- apply(ML.Est.DE.Mat,2,median)
        ML.Est.IE <- apply(ML.Est.IE.Mat,2,median)
        ML.SE.DE <- sqrt(apply((ML.Est.DE.Mat - matrix(ML.Est.DE,100,4,byrow=T))^2 + ML.SE.DE.Mat^2,2,median))
        ML.SE.IE <- sqrt(apply((ML.Est.IE.Mat - matrix(ML.Est.IE,100,4,byrow=T))^2 + ML.SE.IE.Mat^2,2,median))
        
        ################################
        # Summary
        ################################
        
        Eff.IPW.DE[index1,index2] <- (IPW.Est.DE)[4]
        Eff.IPW.IE[index1,index2] <- (IPW.Est.IE)[4]
        Eff.ML.DE[index1,index2] <- (ML.Est.DE)[4]
        Eff.ML.IE[index1,index2] <- (ML.Est.IE)[4]
        
        SE.IPW.DE[index1,index2] <- (IPW.SE.DE)[4]
        SE.IPW.IE[index1,index2] <- (IPW.SE.IE)[4]
        SE.ML.DE[index1,index2] <- (ML.SE.DE)[4]
        SE.ML.IE[index1,index2] <- (ML.SE.IE)[4]
        
        Wald.IPW.DE[index1,index2] <- (IPW.Est.DE/IPW.SE.DE)[4]
        Wald.IPW.IE[index1,index2] <- (IPW.Est.IE/IPW.SE.IE)[4]
        Wald.ML.DE[index1,index2] <- (ML.Est.DE/ML.SE.DE)[4]
        Wald.ML.IE[index1,index2] <- (ML.Est.IE/ML.SE.IE)[4]
        
    }
    
}

Summary <- list()
Summary$Eff.IPW.DE <- Eff.IPW.DE
Summary$Eff.IPW.IE <- Eff.IPW.IE
Summary$Eff.ML.DE  <- Eff.ML.DE
Summary$Eff.ML.IE  <- Eff.ML.IE

Summary$SE.IPW.DE <- SE.IPW.DE
Summary$SE.IPW.IE <- SE.IPW.IE
Summary$SE.ML.DE  <- SE.ML.DE
Summary$SE.ML.IE  <- SE.ML.IE

Summary$Wald.IPW.DE <- Wald.IPW.DE
Summary$Wald.IPW.IE <- Wald.IPW.IE
Summary$Wald.ML.DE  <- Wald.ML.DE
Summary$Wald.ML.IE  <- Wald.ML.IE

# Save Data: this will make Summary.RData file in the current folder.

# save(Summary,file="Summary.RData")

#################################################
# Load Summary and Plot
#################################################

load("Summary.RData")

Eff.IPW.DE <- Summary$Eff.IPW.DE
Eff.IPW.IE <- Summary$Eff.IPW.IE
Eff.ML.DE  <- Summary$Eff.ML.DE
Eff.ML.IE  <- Summary$Eff.ML.IE
SE.IPW.DE <- Summary$SE.IPW.DE
SE.IPW.IE <- Summary$SE.IPW.IE
SE.ML.DE  <- Summary$SE.ML.DE
SE.ML.IE  <- Summary$SE.ML.IE
Wald.IPW.DE <- Summary$Wald.IPW.DE
Wald.IPW.IE <- Summary$Wald.IPW.IE
Wald.ML.DE  <- Summary$Wald.ML.DE
Wald.ML.IE  <- Summary$Wald.ML.IE

#################################################
# Figure 1 : Reletive Efficiency
#################################################


Trt.SC <- 6851 ; Total.SC <- 10907
Trt.SUBA <- 1133 ; Total.SUBA <- 2526
g1 <- (0:20)/20
g2 <- (0:20)/20
alpha2.grid <- c(0,0,Trt.SC/Total.SC,Trt.SUBA/Total.SUBA,Trt.SC/Total.SC,Trt.SUBA/Total.SUBA)

layout(matrix(c(1,2,3),1,3,byrow=T), widths=c(3,3,1), heights=c(1))
par(mar=c(3,3,2,0.5),oma=c(2,2,0,0))

range((SE.IPW.DE/SE.ML.DE)^2)
range((SE.IPW.IE/SE.ML.IE)^2)

ZL <- c(65,110)
ZL2 <- seq(65,110,by=2.5)
ZL3 <- c(70,80,90,100,110)

# ZL <- c(9,41)
# ZL2 <- seq(9,41,by=1)
# ZL3 <- c(9,20,30,41)

par(mar=c(3,3,2,0.5))

filled.contour3( g1, g2, (SE.IPW.DE/SE.ML.DE)^2, color = C2, zlim = ZL, levels=ZL2,
                 plot.axes = {axis(1); axis(2); 
                     # contour(g1, g2, (SE.IPW.DE/SE.ML.DE)^2, levels=ZL3,add=TRUE,lwd=1,col=rgb(0,0,0,0.5));
                     points(alpha2.grid[3],alpha2.grid[4],cex=3,lwd=3,col=4,pch=4)})
title(xlab="",line=2.5)
title(ylab="",line=2.75)
title(main=expression("Direct Effect"),cex.main=1.5)


filled.contour3( g1, g2, (SE.IPW.IE/SE.ML.IE)^2, color = C2, zlim = ZL, levels=ZL2,
                 plot.axes = {axis(1); axis(2); 
                     # contour(g1, g2, (SE.IPW.IE/SE.ML.IE)^2, levels=ZL3,add=TRUE,lwd=1,col=rgb(0,0,0,0.5));
                     points(alpha2.grid[3],alpha2.grid[4],cex=3,lwd=3,col=4,pch=4)})
title(xlab="",line=2.5)
title(ylab="",line=2.75)
title(main=expression("Indirect Effect"),cex.main=1.5)



plot.new()
par(mar=c(3,3,2,3))

filled.legend( c(0,1), c(0,1), matrix(c(9,41,0,0),2,2), color = C2, zlim = ZL, levels=ZL2 )


mtext(expression(alpha["SC"]),side=1,line=0.5,outer=TRUE,cex=1)
mtext(expression(alpha["Su"]),side=2,line=0.5,outer=TRUE,cex=1)



#################################################
# Figure 2 : DE/IE Effects
#################################################

Trt.SC <- 6851 ; Total.SC <- 10907
Trt.SUBA <- 1133 ; Total.SUBA <- 2526
g1 <- (0:20)/20
g2 <- (0:20)/20
alpha2.grid <- c(0,0,Trt.SC/Total.SC,Trt.SUBA/Total.SUBA,Trt.SC/Total.SC,Trt.SUBA/Total.SUBA)

layout(matrix(c(1,2,3),1,3,byrow=T), widths=c(3,3,1))
par(mar=c(3,3,2,0.5),oma=c(2,2,0,0))

ZL <- c(-0.010,0.021)*100
ZL2 <- seq(-0.010,0.021,by=0.001)*100
ZL3 <- c(-qnorm(0.975),0,qnorm(0.975))

range(Eff.ML.DE*100)
range(Eff.ML.IE*100)

L1 <- contourLines(g1,g2,Wald.ML.DE,levels=c(qnorm(0.975)))
kk <- 1
L1x <- c(L1[[kk]]$x[1],L1[[kk]]$x,L1[[kk]]$x[length(L1[[kk]]$x)],0)
L1y <- c(L1[[kk]]$y[1],L1[[kk]]$y,L1[[kk]]$y[length(L1[[kk]]$y)],1)

filled.contour3( g1, g2, Eff.ML.DE*100, color = C1, zlim=ZL, levels=ZL2,
                 plot.axes = {axis(1); axis(2); 
                     # contour(g1, g2, Wald.ML.DE, levels=c(-qnorm(0.975),qnorm(0.975)),add=TRUE,lwd=1,col=rgb(0,0,0,0.5),
                     #         labels=c("",""));
                     polygon(L1x,L1y,border="black",col=rgb(0,0,0,0.0),lwd=0.5);
                     text(0.2,0.8,expression("p-value<0.05"));
                     text(0.8,0.2,expression("p-value">="0.05"));
                     points(alpha2.grid[3],alpha2.grid[4],cex=3,lwd=3,col=4,pch=4)})
title(xlab="",line=2.5)
title(ylab="",line=2.75)
title(main=expression("Direct Effect"),cex.main=1.5)

L2 <- contourLines(g1,g2,Wald.ML.IE,levels=c(qnorm(0.975)))
kk <- 1
L2x <- c(L2[[kk]]$x,L2[[kk]]$x[length(L2[[kk]]$x)])
L2y <- c(L2[[kk]]$y,L2[[kk]]$y[1])
# kk <- 1
# L3x <- c(L2[[kk]]$x)
# L3y <- c(L2[[kk]]$y)

filled.contour3( g1, g2, Eff.ML.IE*100, color = C1, zlim=ZL, levels=ZL2,
                 plot.axes = {axis(1); axis(2); 
                     # contour(g1, g2, Wald.ML.IE, levels=c(-qnorm(0.975),qnorm(0.975)),add=TRUE,lwd=1,col=rgb(0,0,0,0.5),
                     #         labels=c("",""));
                     polygon(L2x,L2y,border="black",col=rgb(0,0,0,0.0),lwd=0.5);
                     text(0.2,0.8,expression("p-value">="0.05"));
                     text(1,0.1,expression("p-value<0.05"),pos=2);
                     points(alpha2.grid[3],alpha2.grid[4],cex=3,lwd=3,col=4,pch=4)})

title(xlab="",line=2.5)
title(ylab="",line=2.75)
title(main=expression("Indirect Effect"),cex.main=1.5)


plot.new()
par(mar=c(3,3,2,3))

filled.legend( c(0,1), c(0,1), matrix(c(-3.5,2.5,0,0),2,2), color = C11, zlim=ZL, levels=ZL2 )

mtext(expression(alpha["SC"]),side=1,line=0.5,outer=TRUE,cex=1)
mtext(expression(alpha["Su"]),side=2,line=0.5,outer=TRUE,cex=1)
