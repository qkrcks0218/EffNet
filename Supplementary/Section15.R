############################################
# Last update : October 02 2021
# Code for the simulation in Section 1.5 in the supplementary material
############################################

############################################
# Packages & Source File
############################################

library(pracma)
source("Function15.R")

##########################
# Parameter Setup
##########################

parameterM <- cbind( rep(c(3,5,7)*0.1,each=99) , rep((1:99)*0.01,3) )   # 1st column: p_A , 2nd column: alpha


# Note: Recommend to use parallel computing 
# by splitting the range of BATCH

for(BATCH in 1:297){
    N <- 10000                      # Number of Clusters
    SigmaX <- 1                     # Variance of Covariates, X ~ N(0,SigmaX^2)
    pA <- parameterM[BATCH,1]       # Randomized Experiment, P(A=1) = pA
    teff <- 3                       # ATE
    beff <- c(1,2,0.5)              # Coefficient of individual outcome regression; g_j^*(a,X_i) = teff*a + beff[1] + beff[2]*X_{ij} + beff[3]*X_{i(-j)}
    Sigma <- c(1,1)                 # Conditional variance of Y_{ij} given a,x
    alpha <- parameterM[BATCH,2]    # Policy parameter
    alpha2 <- 0.9                   # Counterfactual parameter for indirect effect
    iter.max <- 1000                # Maximum iteration
    
    DEestimates <- rep(0,iter.max)  # DE estimates
    IEestimates <- rep(0,iter.max)  # IE estimates
    ATEestimates <- rep(0,iter.max) # ATE estimates
    
    for(iter in 1:iter.max){
        
        ##########################
        # DGP
        ##########################
        
        Xmat <- genX(N,SigmaX)                                          # X matrix (row: each cluster, column: each unit)
        Amat <- genA(Xmat,pA)                                           # A matrix (row: each cluster, column: each unit) + Cluster Propensity Score matrix (00,01,10,11,realized)
        Ymat <- genY(Amat,Xmat,teff,beff,Sigma)                         # Y matrix (row: each cluster, column: each unit)
        
        C <- c(1:N,1:N)                                                 # Cluster variable
        X <- cbind( c(Xmat[,1],Xmat[,2]) , c(Xmat[,2],Xmat[,1]) )       # Making "long format" X matrix (row: each unit)
        A <- cbind( c(Amat[,2],Amat[,3]) , c(Amat[,3],Amat[,2]) )       # Making "long format" A matrix (row: each unit)
        psO <- c(Amat[,8],Amat[,8])                                     # Making "long format" true propensity score (row: each unit)
        Y <- cbind( c(Ymat[,1],Ymat[,2]) )                              # Making "long format" Y vector (row: each unit)
        
        Data <- cbind(Y,A,X,C,psO)                                      # Making data matrix to use in the estimation
        Data <- Data[order(C),]
        Data <- data.frame(Data)
        colnames(Data) <- c("Y","Aj","Ap","Xj","Xp","C","psO")
        
        ##########################
        # DE/IE/ATE estimation
        ##########################
        
        DEestimates[iter] <- DEest(Data,Ymat,Amat,Xmat,alpha,pA=0,peer=1)               # DE estimation
        IEestimates[iter] <- IEest(Data,Ymat,Amat,Xmat,alpha,alpha2,pA=0,peer=1)        # IE estimation
        ATEestimates[iter] <- ATEest(Data,Ymat,Amat,Xmat,pA=0,peer=1)                   # ATE estimation
        
        print(iter)
        
    }
    
    #############################################
    # Result : (iter.max+1) x 3 matrix
    # 1st column: DE
    # 2nd column: ATE
    # 3rd column: IE
    # First iter.max rows: DE/ATE/IE estimates
    # Last row: SEB for DE/ATE/IE 
    #############################################
    
    Result <- cbind( c( DEestimates, SEB.DE(Sigma,SigmaX,pA,alpha) ) , c( ATEestimates, SEB.TRUE.ATE(Sigma,SigmaX,pA) ) , c(IEestimates,0) )        
    write.csv(Result,sprintf("RESULT_pA%03d_al%03d.csv",round(100*pA),round(100*alpha)),row.names=FALSE)
    
}

##########################
# Drawing Plots
##########################

pAlist <- c(0.3,0.5,0.7)            # p_A values
DE.Bound <- matrix(0,99,3)          # SEB for DE
IE.Bound <- matrix(0,99,3)          # SEB for IE

for(ii in 1:99){
    for(jj in 1:3){
        DE.Bound[ii,jj] <- SEB.DE(Sigma,SigmaX,pAlist[jj],ii/100)
        IE.Bound[ii,jj] <- SEB.IE(Sigma,SigmaX,pAlist[jj],ii/100,alpha2)
    }
}

result <- list()                    # Reading all 3*99=297 files
for(pA in 1:3){
    result[[pA]] <- list()
    for(alpha in 1:99){
        result[[pA]][[alpha]] <- read.csv(sprintf("Results/RESULT_pA%03d_al%03d.csv",100*pAlist[pA],alpha))     # Load the Provided Datasets
        print(100*(pA-1)+alpha)
    }
}

TheoBound <- list()                 # Theoretical SEB
EmpBound <- list()                  # Empirical SEB
BiasBound <- list()                 # Empirical Bias

for(pA in 1:3){
    TheoBound[[pA]] <- matrix(0,99,3)
    EmpBound[[pA]] <- matrix(0,99,3)
    BiasBound[[pA]] <- matrix(0,99,3)
    for(alpha in 1:99){
        TheoBound[[pA]][alpha,1] <- result[[pA]][[alpha]][1001,1]/10000
        TheoBound[[pA]][alpha,2] <- result[[pA]][[alpha]][1001,2]/10000
        TheoBound[[pA]][alpha,3] <- result[[pA]][[alpha]][1001,3]/10000
        
        EmpBound[[pA]][alpha,1] <- var(result[[pA]][[alpha]][1:1000,1])
        EmpBound[[pA]][alpha,2] <- var(result[[pA]][[alpha]][1:1000,2])
        EmpBound[[pA]][alpha,3] <- var(result[[pA]][[alpha]][1:1000,3])
        
        BiasBound[[pA]][alpha,1] <- mean(result[[pA]][[alpha]][1:1000,1])-3
        BiasBound[[pA]][alpha,2] <- mean(result[[pA]][[alpha]][1:1000,2])-3
        BiasBound[[pA]][alpha,3] <- mean(result[[pA]][[alpha]][1:1000,3])
    }
}


alpha.grid <- (1:99)*0.01                   # policy parameter grid
yl <- c(-25,25)/10000                       # bias plot range

##########################
# Figure
##########################

mm <- matrix(c(1,2,3,rep(c(4,5,6),5),rep(c(7,8,9),5)),11,3,byrow=T)
layout(mm)

par(mar=c(0,3.5,0,0.5),oma=c(0.5,0,0.1,0.1))

for(pA in 1:3){
    plot.new()
    text(0.5,0.5,bquote(p[A]*"*"*"="*.(pAlist[pA])),cex=1.4)
}

par(mar=c(3,3.5,0.5,0.5))

for(pA in 1:3){
    plot(alpha.grid,BiasBound[[pA]][,1],type='l',col=2, lty=1, lwd=1, ylim=yl,ylab="",yaxt='n')
    if(pA==1){
        title(ylab=expression("Bias"), line=2.25, cex.lab=1.2)
    }
    axis(2,at=c(-0.0025,0,0.0025),label=c("-2.5e-3","0","2.5e-3"))
    abline(h=0,col=1,lty=1,lwd=1.5)
}

yl <- c(0,8)/10000

for(pA in 1:3){
    plot(alpha.grid,TheoBound[[pA]][,1],type='l',lty=1, lwd=1.5, ylim=yl ,ylab="",yaxt='n')
    par(new=T)
    plot(alpha.grid,TheoBound[[pA]][,2],type='l',lty=2, lwd=1.5, col=4, ylim=yl ,ylab="",xaxt='n',yaxt='n')
    par(new=T)
    plot(alpha.grid,EmpBound[[pA]][,1],type='l', lty=1, lwd=1, col=2, ylim=yl ,ylab="",xaxt='n',yaxt='n')
    if(pA==1){
        title(ylab=expression("Variance"), line=2.25, cex.lab=1.2)
    }
    axis(2,at=c(0,0.0004,0.0008),labels=c("0","4e-4","8e-4"))
    points(pAlist[pA],TheoBound[[pA]][1,2]-0.4/10000, pch=17, cex=2)
}
mtext(expression(alpha), side=1, line=-0.5, outer=TRUE, cex=1)



