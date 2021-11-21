############################################
# Last update : October 02 2021
# Code for the data analysis in Section 5
############################################

PI <- function(a,j,alpha){                  # Bernoulli treatment assignment under alpha-strategy
    if(length(a)>1){
        aeij <- a[-j]
        return(prod(alpha^(aeij)*(1-alpha)^(1-aeij)))
    } else {
        return(1)
    }
}

makingAmat <- function(Mk){                 # Returning all possible Mk-dimensional 0/1 Vectors
    if(Mk==1){
        Amat <- matrix(c(0,1),2,1)
    } else if(Mk==2){
        Amat <- cbind(c(0,0,1,1),c(0,1,0,1))
    } else if(Mk==3){
        Amat <- cbind(c(0,0,0,0,1,1,1,1),c(0,0,1,1,0,0,1,1),c(0,1,0,1,0,1,0,1))
    } else if(Mk==4){
        Amat <- cbind(c(rep(0,8),rep(1,8)),
                      c(rep(0,4),rep(1,4),rep(0,4),rep(1,4)),
                      rep(c(0,0,1,1),4),
                      rep(c(0,1),8))
    } else {
        Amat <- cbind(c(rep(0,16),rep(1,16)),
                      c(rep(0,8),rep(1,8),rep(0,8),rep(1,8)),
                      rep(c(0,0,0,0,1,1,1,1),4),
                      rep(c(0,0,1,1),8),
                      rep(c(0,1),16))
    }
    return(Amat)
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

MySL <- function( Data, locY, locX, Ydist=gaussian(), SL.list=c(1:11), MTRY=c(2,4,6,8), MLPL=c(2,4,6,8), MLPdecay=c(10^(-4),10^(-5)), NMN=c(20), obsWeights=NULL ){
    
    ## Poisson c(2,4,5)
    if(Ydist$family=="poisson"){
        SL.list <- intersect(c(1,2,4,5),SL.list)
    }
    
    Learners <- list()
    
    ##############
    # Caret Based
    ##############
    
    SL.caret.SLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      preProc =  c('center', 'scale', 'pca'),
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      learnFuncParams=decay,
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    SL.caret.MLP <- function (Y, X, newX, family, obsWeights, method = "mlpML", 
                              L1,L2,L3,decay,
                              trControl = caret::trainControl(method = "none"), 
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...) 
    {
        if (family$family == "gaussian") {
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      learnFuncParams=decay,
                                      hiddenActFunc = "Act_Logistic",linOut=TRUE,
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights, 
                                      metric = metric, method = method, 
                                      tuneGrid = expand.grid(layer1=L1,layer2=L2,layer3=L3), 
                                      hiddenActFunc = "Act_Identity",
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    
    SL.caret.gbm <- function (Y, X, newX, family, obsWeights, method = "gbm",
                              ntree,intdepth,sh,nmn,
                              trControl = caret::trainControl(method = "none"),
                              metric = ifelse(family$family == "gaussian", "RMSE", "Accuracy"), ...)
    {
        if (family$family == "gaussian") {
            
            fit.train <- caret::train(x = X, y = Y, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "raw")
        }
        if (family$family == "binomial") {
            Y.f <- as.factor(Y)
            levels(Y.f) <- c("A0", "A1")
            fit.train <- caret::train(x = X, y = Y.f, weights = obsWeights,
                                      metric = metric, method = method,
                                      tuneGrid = expand.grid(n.trees=ntree,interaction.depth=intdepth,shrinkage=sh,n.minobsinnode=nmn),
                                      preProc =  c('center', 'scale', 'pca'),
                                      trControl = trControl)
            pred <- predict(fit.train, newdata = newX, type = "prob")[,2]
        }
        fit <- list(object = fit.train)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.caret")
        return(out)
    }
    
    
    #############################################
    # SL-based
    #############################################
    
    SL.new.earth <- function (Y, X, newX, family, obsWeights, id, degree = 2, penalty = 3, 
                              nk = max(21, 2 * ncol(X) + 1), pmethod = "backward", 
                              nfold = 0, ncross = 1, minspan = 0, endspan = 0, ...)
    {
        if (family$family == "gaussian") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights, 
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan)
        }
        if (family$family == "binomial") {
            fit.earth <- earth::earth(x = X, y = Y, degree = degree, weights=obsWeights,
                                      nk = nk, penalty = penalty, pmethod = pmethod, nfold = nfold, 
                                      ncross = ncross, minspan = minspan, endspan = endspan, 
                                      glm = list(family = binomial))
        }
        pred <- predict(fit.earth, newdata = newX, type = "response")
        fit <- list(object = fit.earth)
        out <- list(pred = pred, fit = fit)
        class(out$fit) <- c("SL.earth")
        return(out)
    }
    
    SL.new.xgboost <- function (Y, X, newX, family, obsWeights, id, ntrees = 1000, 
                                max_depth = 4, shrinkage = 0.1, minobspernode = 10, params = list(), 
                                nthread = 1, verbose = 0, save_period = NULL, ...) 
    {
        if (packageVersion("xgboost") < 0.6) 
            stop("SL.xgboost requires xgboost version >= 0.6, try help('SL.xgboost') for details")
        if (!is.matrix(X)) {
            X = model.matrix(~. - 1, X)
        }
        xgmat = xgboost::xgb.DMatrix(data = X, label = Y, weight = obsWeights)
        if (family$family == "gaussian") {
            model = xgboost::xgboost(data = xgmat, objective = "reg:linear", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "binomial") {
            model = xgboost::xgboost(data = xgmat, objective = "binary:logistic", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, nthread = nthread, 
                                     params = params, save_period = save_period)
        }
        if (family$family == "multinomial") {
            model = xgboost::xgboost(data = xgmat, objective = "multi:softmax", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, num_class = length(unique(Y)), 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (family$family == "poisson") {
            model = xgboost::xgboost(data = xgmat, objective = "count:poisson", 
                                     nrounds = ntrees, max_depth = max_depth, min_child_weight = minobspernode, 
                                     eta = shrinkage, verbose = verbose, 
                                     nthread = nthread, params = params, save_period = save_period)
        }
        if (!is.matrix(newX)) {
            newX = model.matrix(~. - 1, newX)
        }
        pred = predict(model, newdata = newX)
        fit = list(object = model)
        class(fit) = c("SL.xgboost")
        out = list(pred = pred, fit = fit)
        return(out)
    }
    
    Learners[[1]] <- create.Learner("SL.glm")
    TOTAL.M <- 1
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.glmnet",tune=list(alpha=c(1,0.5,0),useMin=c(TRUE,FALSE))) 
    # Lasso.min , EN.min , Ridge.min , Lasso.1se , EN.1se , Ridge.1se
    TOTAL.M <- TOTAL.M+1           # 2
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.earth",tune=list(degree=c(1,2,3,4,5)))
    # Earth.deg=1 , Earth.deg=2 , Earth.deg=3 , Earth.deg=4 , Earth.deg=5
    TOTAL.M <- TOTAL.M+1           # 3
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.gam",tune=list(deg.gam=c(1,2,3,4,5)))
    # Gam.deg=1 , Gam.deg=2 , Gam.deg=3 , Gam.deg=4 , Gam.deg=5
    TOTAL.M <- TOTAL.M+1           # 4
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.new.xgboost",tune=list(n.trees=c(100,300,500),max_depth=c(1,2,3,4)))
    # xgboost
    TOTAL.M <- TOTAL.M+1           # 5
    
    # Learners[[5]] <- create.Learner("SL.kernelKnn",tune=list(k=c(1,5,10,20)))
    # knn
    
    # Learners[[5]] <- create.Learner("SL.ksvm",tune=list(kernel=c("rbfdot","polydot","tanhdot")))
    # SVM
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.polymars",tune=list(knots=c(2,3,4)))
    # polspline (similar to earth)
    TOTAL.M <- TOTAL.M+1           # 6
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.ranger",tune=list(num.trees=c(500,1000,1500),mtry=MTRY))
    # RF
    TOTAL.M <- TOTAL.M+1           # 7
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.gbm",tune=list(ntree=c(100,300,500),intdepth=c(1,2,3),sh=c(0.1,0.01),nmn=NMN))
    # gbm
    TOTAL.M <- TOTAL.M+1           # 8
    
    Learners[[TOTAL.M+1]] <- create.Learner("SL.caret.SLP", tune = list(L1=MLPL,L2=c(0),L3=c(0),decay=MLPdecay))
    # 1layer
    TOTAL.M <- TOTAL.M+1           # 9
    
    # Learners[[12]] <- create.Learner("SL.bartMachine",tune=list(num_trees=c(50,100,150),mem_cache_for_speed=FALSE))
    # bart
    
    BaseLearner <- Learners[[ SL.list[1] ]]$names
    
    if( length(SL.list)>1 ){
        for( METHOD in 2:length(SL.list) ){
            BaseLearner <- c(BaseLearner, Learners[[ SL.list[METHOD] ]]$names)
        }
    }
    
    if( length(locX)==1 ){
        dX <- data.frame(matrix( Data[,locX], dim(Data)[1], 1))
        colnames(dX) <- colnames(Data)[locX]
    } else {
        dX <- Data[,locX]
    }
    
    capture.output( Fitted.SL <- SuperLearner(Y=Data[,locY],X=dX,family=Ydist,
                                              SL.library=BaseLearner,cvControl = list(V = 5), obsWeights=obsWeights) , file=NULL )
    
    return(Fitted.SL)
}



filled.contour3 <-
    function (x = seq(0, 1, length.out = nrow(z)),
              y = seq(0, 1, length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes,mar, ...) 
    {
        # modification by Ian Taylor of the filled.contour function
        # to remove the key and facilitate overplotting with contour()
        # further modified by Carey McGilliard and Bridget Ferris
        # to allow multiple plots on one page
        
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                }
                else {
                    z <- x
                    x <- seq.int(0, 1, length.out = nrow(z))
                }
            }
            else stop("no 'z' matrix specified")
        }
        else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
        # mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        # on.exit(par(par.orig))
        # w <- (3 + mar.orig[2]) * par("csi") * 2.54
        # par(las = las)
        # mar <- mar.orig
        plot.new()
        # par(mar=mar)
        plot.window(xlim, ylim, "", xaxs = xaxs, yaxs = yaxs, asp = asp)
        if (!is.matrix(z) || nrow(z) <= 1 || ncol(z) <= 1) 
            stop("no proper 'z' matrix specified")
        if (!is.double(z)) 
            storage.mode(z) <- "double"
        .filled.contour(as.double(x), as.double(y), z, as.double(levels), 
                        col = col)
        if (missing(plot.axes)) {
            if (axes) {
                title(main = "", xlab = "", ylab = "")
                Axis(x, side = 1)
                Axis(y, side = 2)
            }
        }
        else plot.axes
        if (frame.plot) 
            box()
        if (missing(plot.title)) 
            title(...)
        else plot.title
        invisible()
    }


filled.legend <-
    function (x = seq(0, 1, length.out = nrow(z)), y = seq(0, 1, 
                                                           length.out = ncol(z)), z, xlim = range(x, finite = TRUE), 
              ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), 
              levels = pretty(zlim, nlevels), nlevels = 20, color.palette = cm.colors, 
              col = color.palette(length(levels) - 1), plot.title, plot.axes, 
              key.title, key.axes, asp = NA, xaxs = "i", yaxs = "i", las = 1, 
              axes = TRUE, frame.plot = axes, ...) 
    {
        # modification of filled.contour by Carey McGilliard and Bridget Ferris
        # designed to just plot the legend
        if (missing(z)) {
            if (!missing(x)) {
                if (is.list(x)) {
                    z <- x$z
                    y <- x$y
                    x <- x$x
                }
                else {
                    z <- x
                    x <- seq.int(0, 1, length.out = nrow(z))
                }
            }
            else stop("no 'z' matrix specified")
        }
        else if (is.list(x)) {
            y <- x$y
            x <- x$x
        }
        if (any(diff(x) <= 0) || any(diff(y) <= 0)) 
            stop("increasing 'x' and 'y' values expected")
        #  mar.orig <- (par.orig <- par(c("mar", "las", "mfrow")))$mar
        #  on.exit(par(par.orig))
        #  w <- (3 + mar.orig[2L]) * par("csi") * 2.54
        #layout(matrix(c(2, 1), ncol = 2L), widths = c(1, lcm(w)))
        #  par(las = las)
        #  mar <- mar.orig
        #  mar[4L] <- mar[2L]
        #  mar[2L] <- 1
        #  par(mar = mar)
        # plot.new()
        plot.window(xlim = c(0, 1), ylim = range(levels), xaxs = "i", 
                    yaxs = "i")
        rect(0, levels[-length(levels)], 1, levels[-1L], col = col)
        if (missing(key.axes)) {
            if (axes) 
                axis(4)
        }
        else key.axes
        box()
    }



C1 <- function(n,alpha=1){
    l1 <- 11
    l2 <-  n-11+1
    cols <- c(if (l1 > 0) {hsv(h = 180/360, s = seq.int(0.5,if (FALSE){ 0.5/k } else { 0 }, length.out = l1), v = 1, alpha = alpha)} , 
              if (l2 > 1) {hsv(h = 0, s = seq.int(0, 0.75, length.out = l2)[-1L], v = 1, alpha = alpha)} )
    cols
    
}
C11 <- function(n,alpha=1){
    l1 <- 11
    l2 <- n-11+1
    cols <- c(if (l1 > 0) {hsv(h = 180/360, s = seq.int(0.5,if (FALSE){ 0.5/k } else { 0 }, length.out = l1), v = 1, alpha = alpha)} , 
              if (l2 > 1) {hsv(h = 0, s = seq.int(0, 0.75, length.out = l2)[-1L], v = 1, alpha = alpha)} )
    cols
    
}
C2 <- function(n,alpha=1){
    cols <- hsv(h = 240/360, s = seq.int(0.05,0.7, length.out = n), v = 1, alpha = alpha)
    cols
    
}


