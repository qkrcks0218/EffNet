############################################
# Last update : October 02 2021
# Code for the simulation in Section 1.6 in the supplementary material
############################################

############################################
# Functions
############################################

logi <- function(x){                        # Logistic function
    return(1/(1+exp(-x)))
}

PI <- function(a,j,alpha){                  # Bernoulli treatment assignment under alpha-strategy
    aeij <- a[-j]
    return(prod(alpha^(aeij)*(1-alpha)^(1-aeij)))
}

makingAmat <- function(Mk){                 # Returning all possible Mk-dimensional 0/1 Vectors
    if(Mk==2){
        Amat <- cbind(c(0,0,1,1),c(0,1,0,1))
    } else if(Mk==3){
        Amat <- cbind(c(0,0,0,0,1,1,1,1),c(0,0,1,1,0,0,1,1),c(0,1,0,1,0,1,0,1))
    } else {
        Amat <- cbind(c(rep(0,8),rep(1,8)),
                      c(rep(0,4),rep(1,4),rep(0,4),rep(1,4)),
                      rep(c(0,0,1,1),4),
                      rep(c(0,1),8))
    }
    return(Amat)
}

FindEffect <- function(betaXg,Mk,alpha,alpha1,alpha2,truemodel=TRUE,Xmean=c(0.5,0),Cmean=c(0)){    # Calculate True DE/IE effect
    potY <- matrix(0,2^Mk,Mk)
    Amat <- makingAmat(Mk)
    for(ii in 1:(2^Mk)){
        for(jj in 1:Mk){
            if(truemodel==TRUE){
                potY[ii,jj] <- betaXg[1] + Amat[ii,jj]*betaXg[2] + sum(Amat[ii,-jj])*betaXg[3] + 
                    Amat[ii,jj]*Cmean[1]*betaXg[4] + sum(Amat[ii,-jj])*Cmean[1]*betaXg[5] +
                    betaXg[6]*Xmean[1] + betaXg[7]*Xmean[2] + betaXg[8]*Xmean[1]*(Mk-1) + betaXg[9]*Xmean[2]*(Mk-1) + betaXg[10]*Cmean[1]
            } else {
                potY[ii,jj] <- betaXg[1] + Amat[ii,jj]*betaXg[2] + sum(Amat[ii,-jj])*betaXg[3] + sum(betaXg[-(1:3)]*Xmean)
            }
        }
    }
    output <- list()
    output$DE <- sum(DEIEmat(Mk,alpha,alpha1,alpha2)$DE*potY)
    output$IE <- sum(DEIEmat(Mk,alpha,alpha1,alpha2)$IE*potY)
    return(output)
}

DEIEmat <- function(Mk,alpha,alpha1,alpha2){    # Returning all possible weight vectors under policy parameters
    Amat <- makingAmat(Mk)
    wDEmat <- matrix(0,2^Mk,Mk)
    wIEmat <- matrix(0,2^Mk,Mk)
    for(ii in 1:(2^Mk)){
        wDEmat[ii,] <- DEvec(Amat[ii,],alpha)
        wIEmat[ii,] <- IEvec(Amat[ii,],alpha1,alpha2)
    }
    output <- list()
    output$DE <- wDEmat
    output$IE <- wIEmat
    return(output)
}

DEvec <- function(a,alpha){                 # Returning weight vector of DE under a and policy parameter
    Mk <- length(a)
    DEvec <- rep(0,Mk)
    for(jj in 1:Mk){
        DEvec[jj] <- as.numeric(a[jj]==1)*PI(a,jj,alpha) - as.numeric(a[jj]==0)*PI(a,jj,alpha)
    }
    return(DEvec/Mk)
}

IEvec <- function(a,alpha1,alpha2){         # Returning weight vector of IE under a and policy parameters
    Mk <- length(a)
    IEvec <- rep(0,Mk)
    for(jj in 1:Mk){
        IEvec[jj] <- as.numeric(a[jj]==0)*PI(a,jj,alpha1) - as.numeric(a[jj]==0)*PI(a,jj,alpha2)
    }
    return(IEvec/Mk)
}

GenerateX <- function(Nk,Mk,SigmaX,rhoX,Cl.X){   # DGP of X
    Sigma <- diag(rep(SigmaX-rhoX,Mk)) + matrix(rhoX,Mk,Mk)
    return( t(t(chol(Sigma))%*%matrix(rnorm(Nk*Mk),Mk,Nk)) + 0.2*matrix(rep(Cl.X,Mk),Nk,Mk) )
}

ReformX <- function(Xmat,pnumber=2,cpnumber=1){        # Making X matrix where each row constains all covariates in a cluster
    Mk <- dim(Xmat)[2]/(pnumber+cpnumber)
    Nk <- dim(Xmat)[1]
    returnX <- matrix(0,Nk*Mk,pnumber*2+cpnumber)       # 2: ego / others
    for(ii in 1:Nk){
        for(jj in 1:Mk){
            returnX[Mk*(ii-1)+jj,1:pnumber] <- Xmat[ii,(1:pnumber)+(jj-1)*(pnumber+cpnumber)]
            returnX[Mk*(ii-1)+jj,2*pnumber+1:cpnumber] <- Xmat[ii,(pnumber+1:cpnumber)+(jj-1)*(pnumber+cpnumber)]
            exclude <- c((1:pnumber)+(jj-1)*(pnumber+cpnumber), 
                         as.vector( t(matrix( (pnumber+cpnumber)*(1:Mk)-cpnumber+1, Mk, cpnumber ) + matrix(rep(1:cpnumber-1,Mk), Mk, cpnumber, byrow=T)) ) )
            if(Mk==2){
                returnX[Mk*(ii-1)+jj,pnumber+(1:pnumber)] <- Xmat[ii,-exclude]
            } else {
                returnX[Mk*(ii-1)+jj,pnumber+(1:pnumber)] <- apply(matrix(Xmat[ii,-exclude],Mk-1,pnumber,byrow=T),2,sum)
            }
            
        }
    }
    return(returnX)
}

GenerateA <- function(LongXmat,Mk,betaXe,lambda){       # DGP of A
    Nk <- dim(LongXmat)[1]/Mk
    rand.eff <- rep(rnorm(Nk)*sqrt(lambda),each=Mk)
    Avec <- rbinom(dim(LongXmat)[1],1,logi(cbind(1,LongXmat)%*%betaXe+rand.eff))
    AvecO <- rep(0,length(Avec))
    for(ii in 1:Nk){
        for(jj in 1:Mk){
            location <- (ii-1)*Mk + (1:Mk)[-jj]
            AvecO[(ii-1)*Mk+jj] <- sum(Avec[location])
        }
    }
    Avec2 <- matrix(Avec,Nk,Mk,byrow=T)
    result <- list()
    result$ALong <- cbind(Avec,AvecO)
    result$AShort <- Avec2
    return(result)
}

GenerateY <- function(LongAmat,LongAXmat,LongXmat,Mk,betaXg,rho,eta){     # DGP of Y
    Nk <- dim(LongXmat)[1]/Mk
    rand.eff <- rep(rnorm(Nk)*sqrt(rho),each=Mk)
    error <- rnorm(Nk*Mk)*sqrt(eta)
    Yvec <- cbind(1,LongAmat,LongAXmat,LongXmat)%*%betaXg + rand.eff + error
    Ymeanvec <- cbind(1,LongAmat,LongAXmat,LongXmat)%*%betaXg
    YvecO <- rep(0,length(Yvec))
    YmeanvecO <- rep(0,length(Ymeanvec))
    for(ii in 1:Nk){
        for(jj in 1:Mk){
            location <- (ii-1)*Mk + (1:Mk)[-jj]
            YvecO[(ii-1)*Mk+jj] <- sum(Yvec[location])
            YmeanvecO[(ii-1)*Mk+jj] <- sum(Ymeanvec[location])
        }
    }
    Yvec2 <- matrix(Yvec,Nk,Mk,byrow=T)
    Ymeanvec2 <- matrix(Ymeanvec,Nk,Mk,byrow=T)
    result <- list()
    result$YLong <- cbind(Yvec,YvecO)
    result$YShort <- Yvec2
    result$YmeanLong <- cbind(Ymeanvec,YmeanvecO)
    result$YmeanShort <- Ymeanvec2
    return(result)
}

IndPS <- function(Xvec,betaXe,b){               # Individual propensity score given x, k, and a random effect
    return(logi(sum(c(1,Xvec)%*%betaXe)+b))
}

ProdPS <- function(Amat,LongXmat,Mk,betaXe,b,pnumber){      # Cluster propensity score given x, k, and a random effect
    if(Mk==2){
        value <- IndPS(LongXmat[1,],betaXe,b)^Amat[1] * (1-IndPS(LongXmat[1,],betaXe,b))^(1-Amat[1]) * 
            IndPS(LongXmat[2,],betaXe,b)^Amat[2] * (1-IndPS(LongXmat[2,],betaXe,b))^(1-Amat[2])
    } else if(Mk==3){
        value <- IndPS(LongXmat[1,],betaXe,b)^Amat[1] * (1-IndPS(LongXmat[1,],betaXe,b))^(1-Amat[1]) * 
            IndPS(LongXmat[2,],betaXe,b)^Amat[2] * (1-IndPS(LongXmat[2,],betaXe,b))^(1-Amat[2]) * 
            IndPS(LongXmat[3,],betaXe,b)^Amat[3] * (1-IndPS(LongXmat[3,],betaXe,b))^(1-Amat[3])
    } else {
        value <- IndPS(LongXmat[1,],betaXe,b)^Amat[1] * (1-IndPS(LongXmat[1,],betaXe,b))^(1-Amat[1]) * 
            IndPS(LongXmat[2,],betaXe,b)^Amat[2] * (1-IndPS(LongXmat[2,],betaXe,b))^(1-Amat[2]) * 
            IndPS(LongXmat[3,],betaXe,b)^Amat[3] * (1-IndPS(LongXmat[3,],betaXe,b))^(1-Amat[3]) * 
            IndPS(LongXmat[4,],betaXe,b)^Amat[4] * (1-IndPS(LongXmat[4,],betaXe,b))^(1-Amat[4])
    }
    return(value)
}

truePS <- function(Amat,LongXmat,Mk,betaXe,lambda,pnumber=2){       # Cluster propensity score given x, k (i.e., integrating out a random effect)
    if(lambda>0.001){
        f <- function(x){ ProdPS(Amat,LongXmat,Mk,betaXe,x,pnumber)*dnorm(x,mean=0,sd=sqrt(lambda)) }
        return(integrate(f,lower=-Inf,upper=Inf)$value)
    } else {
        return( ProdPS(Amat,LongXmat,Mk,betaXe,0,pnumber) )
    }
    
}

ResidualWeighting <- function(Residualmat,Amat,PSmat,Mk,alpha,alpha1,alpha2){       # First part of EIF
    temp <- DEIEmat(Mk,alpha,alpha1,alpha2)
    wDEmat <- temp$DE
    wIEmat <- temp$IE
    
    Nk <- dim(Residualmat)[1]
    A.DEmat <- matrix(0,Nk,Mk)
    A.IEmat <- matrix(0,Nk,Mk)
    if(Mk==2){
        aindex <- 2*Amat[,1]+Amat[,2]+1
    } else if(Mk==3) {
        aindex <- 4*Amat[,1]+2*Amat[,2]+Amat[,3]+1
    } else {
        aindex <- 8*Amat[,1]+4*Amat[,2]+2*Amat[,3]+Amat[,4]+1
    }
    
    for(ii in 1:Nk){
        A.DEmat[ii,] <- wDEmat[aindex[ii],]
        A.IEmat[ii,] <- wIEmat[aindex[ii],]
    }
    output <- list()
    output$DE <- apply(A.DEmat*Residualmat,1,sum)/PSmat
    output$IE <- apply(A.IEmat*Residualmat,1,sum)/PSmat
    return(output)
}

ORmean <- function(LongXmat,position.of.interaction,position.of.X,Mk,betaXg,alpha,alpha1,alpha2,pnumber=2,cpnumber=1){         # Second part of EIF
    
    Nk <- dim(LongXmat)[1]/Mk
    DE.ORest <- rep(0,Nk)
    IE.ORest <- rep(0,Nk)
    
    Amat <- makingAmat(Mk)
    
    temp <- DEIEmat(Mk,alpha,alpha1,alpha2)
    wDEmat <- temp$DE
    wIEmat <- temp$IE
    
    if(is.null(position.of.interaction)){
        
        for(member in 1:Nk){
            EstY <- matrix(0,2^Mk,Mk)
            for(ii in 1:(2^Mk)){
                for(jj in 1:Mk){
                    
                    egoA <- Amat[ii,jj] ; otherA <- sum(Amat[ii,-jj])
                    tempX <- (LongXmat[ (member-1)*Mk + jj, position.of.X])
                    
                    EstY[ii,jj] <- betaXg[1] + egoA*betaXg[2] + otherA*betaXg[3] + 
                        sum( tempX*betaXg[3+(1:length(tempX))] )
                }
            }
            DE.ORest[member] <- sum(wDEmat*EstY)
            IE.ORest[member] <- sum(wIEmat*EstY)
        }
        
    } else {
        
        for(member in 1:Nk){
            EstY <- matrix(0,2^Mk,Mk)
            for(ii in 1:(2^Mk)){
                for(jj in 1:Mk){
                    
                    egoA <- Amat[ii,jj] ; otherA <- sum(Amat[ii,-jj])
                    tempAX <- c(LongXmat[ (member-1)*Mk + jj, ]*egoA , LongXmat[ (member-1)*Mk + jj, ]*otherA)[position.of.interaction]
                    tempX <- (LongXmat[ (member-1)*Mk + jj, position.of.X])
                    
                    EstY[ii,jj] <- betaXg[1] + egoA*betaXg[2] + otherA*betaXg[3] + 
                        sum( tempAX*betaXg[3+(1:length(tempAX))] ) +
                        sum( tempX*betaXg[3+length(tempAX)+(1:length(tempX))] )
                }
            }
            DE.ORest[member] <- sum(wDEmat*EstY)
            IE.ORest[member] <- sum(wIEmat*EstY)
        }
        
    }
    
    
    output <- list()
    output$DE <- DE.ORest
    output$IE <- IE.ORest
    
    return(output)
}
