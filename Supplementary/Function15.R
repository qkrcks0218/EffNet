############################################
# Last update : October 02 2021
# Code for the simulation in Section 1.5 in the supplementary material
############################################

############################################
# Functions
############################################

logi <- function(x){                                    # Logistic function
    return(1/(1+exp(-x)))
}

policy <- function(j,avec,alpha){                       # Bernoulli treatment assignment under alpha-strategy
    return(alpha^(avec[-j])*(1-alpha)^(1-avec[-j]))
}

genX <- function(N,SigmaX){                             # DGP of X
    X1 <- rnorm(N)*SigmaX
    X2 <- rnorm(N)*SigmaX
    return( cbind(X1,X2) )
}

genA <- function(Xmat,pA){                              # DGP of A
    N <- dim(Xmat)[1]
    avec1 <- rbinom(N,1,pA)
    avec2 <- rbinom(N,1,pA)
    avec <- 2*avec1 + avec2 + 1
    ps00 <- rep((1-pA)^2,N)
    ps01 <- rep((1-pA)*pA,N)
    ps10 <- rep((1-pA)*pA,N)
    ps11 <- rep(pA^2,N)
    psobs <- pA^(avec1+avec2)*(1-pA)^(2-avec1-avec2)
    
    return( cbind(avec,avec1,avec2,ps00,ps01,ps10,ps11,psobs) )
}

IOR <- function(j,avec,xvec,teff,beff){                 # Individual OR model
    return( avec[j]*teff + sum(cbind(1,xvec[j],xvec[-j])*beff ) )
}

genY <- function(Amat,Xmat,teff,beff,Sigma){            # DGP of Y
    N <- dim(Xmat)[1]
    Y1 <- rep(0,N)
    Y2 <- rep(0,N)
    IOR1 <- rep(0,N)
    IOR2 <- rep(0,N)
    for(ii in 1:N){
        Y1[ii] <- IOR(1,Amat[ii,2:3],Xmat[ii,],teff,beff) + rnorm(1)*Sigma[1]
        Y2[ii] <- IOR(2,Amat[ii,2:3],Xmat[ii,],teff,beff) + rnorm(1)*Sigma[2]
        IOR1[ii] <- IOR(1,Amat[ii,2:3],Xmat[ii,],teff,beff)
        IOR2[ii] <- IOR(2,Amat[ii,2:3],Xmat[ii,],teff,beff)
    }
    return( cbind(Y1,Y2,IOR1,IOR2) )
}

genH <- function(PS,Sigma){                                             # H matrix in Lemma A.4 of supplemenetary material 
    if(length(dim(PS))>0){
        N <- dim(PS)[1]
        ps00 <- PS[,1]
        ps01 <- PS[,2]
        ps10 <- PS[,3]
        ps11 <- PS[,4]
        H00 <- matrix(0,dim(PS)[1],2)
        H01 <- matrix(0,dim(PS)[1],2)
        H10 <- matrix(0,dim(PS)[1],2)
        H11 <- matrix(0,dim(PS)[1],2)
        for(ii in 1:N){
            Cmat <- matrix(c(0,0,0,0,ps10[ii],0,ps11[ii],0,
                             ps00[ii],0,ps01[ii],0,0,0,0,0,
                             0,0,0,ps01[ii],0,0,0,ps11[ii],
                             0,ps00[ii],0,0,0,ps10[ii],0,0,
                             0,0,0,0,-Sigma[1],0,Sigma[1],0,
                             Sigma[1],0,-Sigma[1],0,0,0,0,0,
                             0,0,0,-Sigma[2],0,0,0,Sigma[2],
                             0,Sigma[2],0,0,0,-Sigma[2],0,0),8,8,byrow=TRUE)
            dmat <- matrix(c(1,-1,1,-1,0,0,0,0),8,1)
            Hmat <- solve(Cmat)%*%dmat
            H00[ii,] <- Hmat[1:2]
            H01[ii,] <- Hmat[3:4]
            H10[ii,] <- Hmat[5:6]
            H11[ii,] <- Hmat[7:8]
        }
        return(cbind(H00,H01,H10,H11))
    } else {
        ps00 <- PS[1]
        ps01 <- PS[2]
        ps10 <- PS[3]
        ps11 <- PS[4]
        Cmat <- matrix(c(0,0,0,0,ps10,0,ps11,0,
                         ps00,0,ps01,0,0,0,0,0,
                         0,0,0,ps01,0,0,0,ps11,
                         0,ps00,0,0,0,ps10,0,0,
                         0,0,0,0,-Sigma[1],0,Sigma[1],0,
                         Sigma[1],0,-Sigma[1],0,0,0,0,0,
                         0,0,0,-Sigma[2],0,0,0,Sigma[2],
                         0,Sigma[2],0,0,0,-Sigma[2],0,0),8,8,byrow=TRUE)
        dmat <- matrix(c(1,-1,1,-1,0,0,0,0)/2,8,1)
        Hmat <- solve(Cmat)%*%dmat
        return(Hmat)
    }
}


quadratic <- function(Sigma,pA){                                    # Coefficient of EIF
    psvec <- c((1-pA)^2,pA*(1-pA),pA*(1-pA),pA^2)
    Hmat <- genH(psvec,Sigma)
    return(psvec[1]*t(Hmat[1:2])%*%diag(Sigma)%*%(Hmat[1:2]) +
               psvec[2]*t(Hmat[3:4])%*%diag(Sigma)%*%(Hmat[3:4]) +
               psvec[3]*t(Hmat[5:6])%*%diag(Sigma)%*%(Hmat[5:6]) +
               psvec[4]*t(Hmat[7:8])%*%diag(Sigma)%*%(Hmat[7:8]) )
}

SEB.TRUE.ATE <- function(Sigma,SigmaX,pA){                          # SEB for ATE calculation
    return( 1/4* (Sigma[1]/pA + Sigma[1]/(1-pA) + Sigma[2]/pA + Sigma[2]/(1-pA)) )
}

SEB.DE <- function(Sigma,SigmaX,pA,alpha){                          # SEB for DE calculation
    Q <- function(x,y){
        psvec <- c((1-pA)^2,pA*(1-pA),pA*(1-pA),pA^2)
        return( ( (policy(1,c(0,0),alpha)^2*Sigma[1] + policy(2,c(0,0),alpha)^2*Sigma[2])/psvec[1] +
                      (policy(1,c(0,1),alpha)^2*Sigma[1] + policy(2,c(0,1),alpha)^2*Sigma[2])/psvec[2] +
                      (policy(1,c(1,0),alpha)^2*Sigma[1] + policy(2,c(1,0),alpha)^2*Sigma[2])/psvec[3] +
                      (policy(1,c(1,1),alpha)^2*Sigma[1] + policy(2,c(1,1),alpha)^2*Sigma[2])/psvec[4] )/4*
                    dnorm(x,0,SigmaX)* dnorm(y,0,SigmaX) )
    }
    return(integral2(Q,xmin=-100,xmax=100,ymin=-100,ymax=100)$Q)
}

SEB.IE <- function(Sigma,SigmaX,pA,alpha1,alpha2){                  # SEB for IE calculation
    Q <- function(x,y){
        psvec <- c((1-pA)^2,pA*(1-pA),pA*(1-pA),pA^2)
        return(  ((Sigma[1]+Sigma[2])/psvec[1]+Sigma[1]/psvec[2]+Sigma[2]/psvec[3])*(alpha1-alpha2)^2/4*dnorm(x,0,SigmaX)* dnorm(y,0,SigmaX) ) 
    }
    return(integral2(Q,xmin=-100,xmax=100,ymin=-100,ymax=100)$Q)
}

DEest <- function(Data,Ymat,Amat,Xmat,alpha,pA=0,peer=1){               # DE estimation
    
    N <- dim(Ymat)[1]
    
    wDE <- matrix(0,2,4)
    wDE[,1] <- c(-policy(1,c(0,0),alpha),-policy(2,c(0,0),alpha))/2
    wDE[,2] <- c(-policy(1,c(0,1),alpha),policy(2,c(0,1),alpha))/2
    wDE[,3] <- c(policy(1,c(1,0),alpha),-policy(2,c(1,0),alpha))/2
    wDE[,4] <- c(policy(1,c(1,1),alpha),policy(2,c(1,1),alpha))/2
    
    
    LMfit <- lm(Y~Aj+Ap+Xj+Xp,data=Data)
    teffe <- LMfit$coefficients[2:3]
    beffe <- LMfit$coefficients[c(1,4,5)]
    Yfit <- matrix(0,N,2)
    Yfit[,1] <- cbind(Amat[,2],Amat[,3])%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit[,2] <- cbind(Amat[,3],Amat[,2])%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    residual <- Ymat[,1:2] - Yfit
    
    Yfit00 <- matrix(0,N,2)
    Yfit00[,1] <- cbind(rep(0,N),rep(0,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit00[,2] <- cbind(rep(0,N),rep(0,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit01 <- matrix(0,N,2)
    Yfit01[,1] <- cbind(rep(0,N),rep(1,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit01[,2] <- cbind(rep(1,N),rep(0,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit10 <- matrix(0,N,2)
    Yfit10[,1] <- cbind(rep(1,N),rep(0,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit10[,2] <- cbind(rep(0,N),rep(1,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit11 <- matrix(0,N,2)
    Yfit11[,1] <- cbind(rep(1,N),rep(1,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit11[,2] <- cbind(rep(1,N),rep(1,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    
    pAhat <- mean(Data$Aj)
    Atotal <- apply(Amat[,2:3],1,sum)
    PShat <- pAhat^(Atotal)*(1-pAhat)^(2-Atotal)
    
    Q1 <- rep(0,N) ; Q2 <- rep(0,N)
    
    for(ii in 1:N){
        Q1[ii] <- sum(wDE[,Amat[ii,1]]*residual[ii,])/PShat[ii]
        Q2[ii] <- sum(wDE[,1]*Yfit00[ii,]) + sum(wDE[,2]*Yfit01[ii,]) + sum(wDE[,3]*Yfit10[ii,]) + sum(wDE[,4]*Yfit11[ii,])
    }
    
    return(sum(Q1+Q2)/N)
}

IEest <- function(Data,Ymat,Amat,Xmat,alpha,alpha2,pA=0,peer=1){              # IE estimation
    
    N <- dim(Ymat)[1]
    
    wIE <- matrix(0,2,4)
    wIE[,1] <- c(policy(1,c(0,0),alpha)-policy(1,c(0,0),alpha2),policy(2,c(0,0),alpha)-policy(2,c(0,0),alpha2))/2
    wIE[,2] <- c(policy(1,c(0,1),alpha)-policy(1,c(0,1),alpha2),0)/2
    wIE[,3] <- c(0,policy(2,c(1,0),alpha)-policy(2,c(1,0),alpha2))/2
    wIE[,4] <- c(0,0)/2
    
    
    LMfit <- lm(Y~Aj+Ap+Xj+Xp,data=Data)
    teffe <- LMfit$coefficients[2:3]
    beffe <- LMfit$coefficients[c(1,4,5)]
    Yfit <- matrix(0,N,2)
    Yfit[,1] <- cbind(Amat[,2],Amat[,3])%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit[,2] <- cbind(Amat[,3],Amat[,2])%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    residual <- Ymat[,1:2] - Yfit
    
    Yfit00 <- matrix(0,N,2)
    Yfit00[,1] <- cbind(rep(0,N),rep(0,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit00[,2] <- cbind(rep(0,N),rep(0,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit01 <- matrix(0,N,2)
    Yfit01[,1] <- cbind(rep(0,N),rep(1,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit01[,2] <- cbind(rep(1,N),rep(0,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit10 <- matrix(0,N,2)
    Yfit10[,1] <- cbind(rep(1,N),rep(0,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit10[,2] <- cbind(rep(0,N),rep(1,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit11 <- matrix(0,N,2)
    Yfit11[,1] <- cbind(rep(1,N),rep(1,N))%*%teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit11[,2] <- cbind(rep(1,N),rep(1,N))%*%teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    
    
    pAhat <- mean(Data$Aj)
    Atotal <- apply(Amat[,2:3],1,sum)
    PShat <- pAhat^(Atotal)*(1-pAhat)^(2-Atotal)
    
    Q1 <- rep(0,N) ; Q2 <- rep(0,N)
    
    for(ii in 1:N){
        Q1[ii] <- sum(wIE[,Amat[ii,1]]*residual[ii,])/PShat[ii]
        Q2[ii] <- sum(wIE[,1]*Yfit00[ii,]) + sum(wIE[,2]*Yfit01[ii,]) + sum(wIE[,3]*Yfit10[ii,]) + sum(wIE[,4]*Yfit11[ii,])
    }
    
    return(sum(Q1+Q2)/N)
}


ATEest <- function(Data,Ymat,Amat,Xmat,pA=0,peer=1){                    # ATE estimation
    
    N <- dim(Ymat)[1]
    
    
    LMfit_noint <- lm(Y~Aj+Xj+Xp,data=Data)
    teffe <- LMfit_noint$coefficients[2]
    beffe <- LMfit_noint$coefficients[c(1,3,4)]
    
    Yfit <- matrix(0,N,2)
    Yfit[,1] <- Amat[,2]*teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit[,2] <- Amat[,3]*teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    residual <- Ymat[,1:2] - Yfit
    
    Yfit00 <- matrix(0,N,2)
    Yfit00[,1] <- rep(0,N)*teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit00[,2] <- rep(0,N)*teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit01 <- matrix(0,N,2)
    Yfit01[,1] <- rep(0,N)*teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit01[,2] <- rep(1,N)*teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit10 <- matrix(0,N,2)
    Yfit10[,1] <- rep(1,N)*teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit10[,2] <- rep(0,N)*teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    Yfit11 <- matrix(0,N,2)
    Yfit11[,1] <- rep(1,N)*teffe + cbind(1,Xmat[,1],Xmat[,2])%*%beffe
    Yfit11[,2] <- rep(1,N)*teffe + cbind(1,Xmat[,2],Xmat[,1])%*%beffe
    
    
    pAhat <- mean(Data$Aj)
    
    Atotal <- apply(Amat[,2:3],1,sum)
    PShat <- pAhat^(Atotal)*(1-pAhat)^(2-Atotal)
    PShatmat <- cbind( rep((1-pAhat)^2,N) , rep((1-pAhat)*(pAhat),N) , rep((1-pAhat)*(pAhat),N) , rep((pAhat)^2,N) )
    
    Sigmahat <- c(var(residual[,1]),var(residual[,2]))*(N-1)/N
    
    Hmat <- genH(PShatmat,Sigmahat)
    
    Q1 <- rep(0,N) ; Q2 <- rep(0,N)
    
    for(ii in 1:N){
        pos <- c(2*Amat[ii,1]-1,2*Amat[ii,1])
        Q1[ii] <- sum(Hmat[ii,pos]*residual[ii,])
        Q2[ii] <- sum(Yfit11[ii,] - Yfit00[ii,])/2
    }
    
    return(sum(Q1+Q2)/N)
}

