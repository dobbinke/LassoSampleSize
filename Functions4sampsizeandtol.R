

# functions for sample size and tolerance estimation

get.params4tolerance <- function(
	OptAcc,
	MCruns,
	p1,
	NmbrInfFeatures,
	CovMatStruct,
	Block.size,
	CorrParamVal,
	Dimension,
	trim.prop,
	type.measure.for.cvglmnet,
	number.of.computer.clusters) {

	my.return <- list(
	OptAcc=OptAcc,
	MCruns=MCruns,
	p1=p1,
	NmbrInfFeatures=NmbrInfFeatures,
	CovMatStruct=CovMatStruct,
	Block.size=Block.size,
	CorrParamVal=CorrParamVal,
	Dimension=Dimension,
	trim.prop=trim.prop,
	type.measure.for.cvglmnet=type.measure.for.cvglmnet,
	number.of.computer.clusters = number.of.computer.clusters		
	)

	my.return;
}


get.covmat <- function(
	NmbrInfFeatures,
	CovMatStruct,
	Block.size,
	CorrParamVal,
	Dimension) {

     condition.check = ((NmbrInfFeatures %% Block.size) == 0)
     if (condition.check == FALSE) {
     	stop("NmbrInfFeatures must be a multiple of Block.size.")
     }

    number.of.blocks <- NmbrInfFeatures/Block.size;
    my.blocks <- rep(1:number.of.blocks,each=Block.size);

	my.return <- diag(rep(1,Dimension))
	if (CovMatStruct == "AR1.block") {
		for (i in 1:NmbrInfFeatures) {
			for (j in 1:NmbrInfFeatures) {
				my.return[i,j] <- (CorrParamVal^abs(i-j))*(my.blocks[i]==my.blocks[j]);
			}
		}
	}

	my.return;
}


get.mulength.TransformedSpace <- function(p1,OptAcc)  {
	Optimal.cutpoint = 0.5 * log( (1-p1)/p1 );
	myfun4root <- function(mulength) {
		OptAcc - p1 *pnorm(mulength-Optimal.cutpoint) -
		(1-p1) * pnorm(mulength+Optimal.cutpoint)
	}
	mymulength <- uniroot(myfun4root,lower=0,upper=10)$root
	mymulength
}


install.packages("matrixcalc")
library(matrixcalc)
get.meanvec.UntransformedSpace <- function(
	NmbrInfFeatures,
	mulength.TransformedSpace,
	covmat.Inverse) {

	my.fun.4.root <- function(xreal) {
		myvec <- xreal * rep(1,NmbrInfFeatures);
		covmat.Inverse.Submatrix <- covmat.Inverse[1:NmbrInfFeatures,1:NmbrInfFeatures]
		myQuad <- sqrt(as.numeric(t(myvec) %*% covmat.Inverse.Submatrix %*% myvec));
		# THE LENGTH OF MU IS HALF THE MAH.D.
		myReturn <- myQuad - (mulength.TransformedSpace) 
		myReturn;
	}
	myroot <- uniroot(my.fun.4.root,lower = 0.001, upper = 10)
	myroot
	meanvec.UntransformedSpace <- myroot$root * rep(1,NmbrInfFeatures)
	the.Dimension <- dim(covmat.Inverse)[1]
	meanvec.UntransformedSpace <- 
	  c(meanvec.UntransformedSpace,
	  	rep(0,the.Dimension-length(meanvec.UntransformedSpace)))
	meanvec.UntransformedSpace
}

get.meanvec.TransformedSpace <- function(
	meanvec.UntransformedSpace,
	covmat.Inverse.Root,
	NmbrInfFeatures) {
	the.Dimension <- dim(covmat.Inverse.Root)[1]
	meanvec.TransformedSpace <- covmat.Inverse.Root %*% meanvec.UntransformedSpace 
 	meanvec.TransformedSpace
}

get.betainfty <- function(meanvec.TransformedSpace,covmat.Inverse)  {
	beta.infinity <- 2 * sqrt(t(meanvec.TransformedSpace) %*% covmat.Inverse %*% meanvec.TransformedSpace)
	beta.infinity
}
install.packages("glmnet")
install.packages("foreach")
install.packages("doParallel")
install.packages("MASS")
library(MASS)
library(glmnet)
library(foreach)
library(doParallel)

EstimateTol <- function(
	current.N,
	MCruns,
    covmat,
    Dimension,
    covmat.Root,
    type.measure,
    number.of.computer.clusters,
    p1,
	meanvec.UntransformedSpace,
	AccuracyORLogistic,
	NmbrInfFeatures,
	OptAcc
	) {

	current.N = current.N
	MCruns = MCruns
    covmat = covmat
    Dimension = Dimension
    covmat.Root = covmat.Root
    type.measure = type.measure
    number.of.computer.clusters = number.of.computer.clusters
    p1 = p1
	meanvec.UntransformedSpace = meanvec.UntransformedSpace
	AccuracyORLogistic = AccuracyORLogistic
	NmbrInfFeatures = NmbrInfFeatures
	OptAcc = OptAcc

	MCxmat <- matrix(NA,nrow=current.N,ncol=Dimension);

    mycluster=makeCluster(number.of.computer.clusters);
    registerDoParallel(mycluster);
    getDoParWorkers();

    myres <- foreach(i=1:MCruns,
    	.combine=rbind,
    	.export=ls()) %dopar% {
    	library(MASS)
    	library(glmnet)

 		MCyresp <- rbinom(current.N,1,p1);
    	MCxmat  <- mvrnorm(n=current.N,mu=meanvec.UntransformedSpace,Sigma=covmat) *2.0*((MCyresp==1)-0.5) 
    	mycvglmnet <- cv.glmnet(x=MCxmat,y=MCyresp,type.measure=type.measure,
            family="binomial",nfolds=10)
        optlambda = mycvglmnet$lambda.min;
        myglmnet = glmnet(x=MCxmat,y=MCyresp,family="binomial",lambda=optlambda)
        myL.UntransformedSpace = as.numeric(coef(myglmnet))
        cutpoint.UntransformedSpace <- myL.UntransformedSpace[1]

        # myL a.k.a. Lhat
        myL.UntransformedSpace = myL.UntransformedSpace[2:(Dimension+1)];
        # n.b. full beta contains intercept, which we need to strip
        # out before standardizing to unit length.

        myL.TransformedSpace = myL.UntransformedSpace %*% covmat.Root;

        myLUnitLength.TransformedSpace <- myL.TransformedSpace/sqrt(sum(myL.TransformedSpace^2))

        return(c(cutpoint.UntransformedSpace,myLUnitLength.TransformedSpace))
    }
    stopCluster(mycluster);

    bigLmatLengths1.TransformedSpace  <- matrix(NA,nrow=MCruns,ncol=Dimension);
    cutpoints.TransformedSpace <- rep(NA,MCruns)

    for (thismc in 1:MCruns) {
        cutpoints.TransformedSpace[thismc] <- myres[thismc,1];
        bigLmatLengths1.TransformedSpace[thismc,] <- myres[thismc,2:(Dimension+1)];
    }

    covmat.Inverse <- matrix.inverse(covmat)
    mulength.TransformedSpace <- get.mulength.TransformedSpace(p1,OptAcc)
    mvu <- get.meanvec.UntransformedSpace(
    		NmbrInfFeatures,
    		mulength.TransformedSpace,
    		covmat.Inverse)
    mvt <- get.meanvec.TransformedSpace(
    	meanvec.UntransformedSpace = mvu,
    	covmat.Inverse.Root,
    	NmbrInfFeatures)
	beta.infinity <- get.betainfty(mvt,covmat.Inverse) 
    if (AccuracyORLogistic=="Accuracy")  {
        theseaccs <- 
          p1 * pnorm(bigLmatLengths1.TransformedSpace %*% mvt - cutpoints.TransformedSpace) +
          (1-p1) * pnorm(bigLmatLengths1.TransformedSpace %*% mvt + cutpoints.TransformedSpace) 
    
        thesetols <- rep(OptAcc,length(theseaccs)) - theseaccs;
    }
    if (AccuracyORLogistic=="Logistic")  {
        myangles <- rep(NA,MCruns);
        costhetas <- rep(NA,MCruns)
        for (i in 1:MCruns) {
            this.beta <- bigLmatLengths1.TransformedSpace[i,]
            this.numerator <- sum(mvt*this.beta);
            this.denominator <- sqrt( sum(mvt*mvt) * sum(this.beta*this.beta))
            myangles[i] <- acos(this.numerator/this.denominator)
            costhetas[i] <- this.numerator/this.denominator
        }
        # Note: cutpoints not used in logistic slope calculation
		beta.infinity <- 2 * sqrt(sum(mvt*mvt)) * 
           sqrt( (1+p1*(1-p1)*sum(mvt*mvt))/(1+p1*(1-p1)*sum(mvt*mvt))) 
        print("Beta infinity is:")
        print(beta.infinity)
        beta.Lhats <- 2 * sqrt(sum(mvt*mvt)) * costhetas * 
           sqrt( (1+p1*(1-p1)*sum(mvt*mvt)*costhetas)/(1+p1*(1-p1)*sum(mvt*mvt)))
        thesetols <- rep(beta.infinity,length(beta.Lhats)) - beta.Lhats;

    }

    summary(thesetols)
    trim.prop = 0;

    returndf <- data.frame(tolest = mean(thesetols,trim=trim.prop), tolestsd = sqrt(var(thesetols)));

    return(returndf);

 }



SampSizeFun <- function(
    TolTarget,
  	LB.N,
  	UB.N,
  	meanvec.UntransformedSpace,
  	MCruns,
  	maxit,
    covmat,
    Dimension,
    covmat.Root,
    type.measure,
    number.of.computer.clusters,
    p1,
    AccuracyORLogistic,
    NmbrInfFeatures,
    OptAcc
	)  {

    meanvec.UntransformedSpace = meanvec.UntransformedSpace;

    # Save the search path N's and tolests 
    myresultsTable <<- data.frame(N.path=rep(NA,maxit),tolest.path=rep(NA,maxit),tolestses=rep(NA,maxit));
    
    myreturn.LB <- EstimateTol(
            current.N = LB.N,
            MCruns = MCruns,
            covmat = covmat,
            Dimension = Dimension,
            covmat.Root = covmat.Root,
            type.measure = type.measure,
            number.of.computer.clusters = number.of.computer.clusters,
            p1 = p1,
            meanvec.UntransformedSpace = meanvec.UntransformedSpace,
            AccuracyORLogistic = AccuracyORLogistic,
            NmbrInfFeatures = NmbrInfFeatures, 
            OptAcc = OptAcc
            )
    myresultsTable$N.path[1] <<- LB.N;
    myresultsTable$tolest.path[1] <<- myreturn.LB$tolest;
    myresultsTable$tolestses[1] <<- myreturn.LB$tolestsd/sqrt(MCruns)
    print(myresultsTable[1,]) 
    
    if (myreturn.LB$tolest < TolTarget) { stop("Smallest sample size achieves tolerance.") }
    
    myreturn.UB <- EstimateTol(
            current.N = UB.N,
            MCruns = MCruns,
            covmat = covmat,
            Dimension = Dimension,
            covmat.Root = covmat.Root,
            type.measure = type.measure,
            number.of.computer.clusters = number.of.computer.clusters,
            p1 = p1,
            meanvec.UntransformedSpace = meanvec.UntransformedSpace,
            AccuracyORLogistic = AccuracyORLogistic,
            NmbrInfFeatures = NmbrInfFeatures, 
            OptAcc = OptAcc
            )
    myresultsTable$N.path[2] <<- UB.N;
    myresultsTable$tolest.path[2] <<- myreturn.UB$tolest;
    myresultsTable$tolestses[2] <<- myreturn.UB$tolestsd/sqrt(MCruns)
    print(myresultsTable[2,]) 
    
    
    if (myreturn.UB$tolest > TolTarget) { stop("Largest sample size cannot achieve tolerance.") }
    
    current.positive.N <- LB.N
    current.negative.N <- UB.N

    if (maxit <3) { stop("Parameter maxit must be at least 3.") }

    last.iter <<- maxit;
    
    for (thisiter in 3:maxit)  {
    
        this.N = round(0.5*(current.positive.N+current.negative.N));
    
        this.object <- EstimateTol(
            current.N = this.N,
            MCruns = MCruns,
            covmat = covmat,
            Dimension = Dimension,
            covmat.Root = covmat.Root,
            type.measure = type.measure,
            number.of.computer.clusters = number.of.computer.clusters,
            p1 = p1,
            meanvec.UntransformedSpace = meanvec.UntransformedSpace,
            AccuracyORLogistic = AccuracyORLogistic,
            NmbrInfFeatures = NmbrInfFeatures, 
            OptAcc = OptAcc
            )
        if ((this.object$tolest-TolTarget) < 0) {
            current.negative.N = this.N
        }
        else {
            current.positive.N = this.N
        }

	    myresultsTable$N.path[thisiter] <<- this.N;
	    myresultsTable$tolest.path[thisiter] <<- this.object$tolest;
	    myresultsTable$tolestses[thisiter] <<- this.object$tolestsd/sqrt(MCruns)
	    print(myresultsTable[thisiter,]) 
    
        if (myresultsTable$N.path[thisiter]==myresultsTable$N.path[thisiter-1]) {
            last.iter <<- thisiter;
            break;
        }
    }
    return(current.positive.N)
}
