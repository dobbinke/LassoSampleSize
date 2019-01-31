




#source("C:/Users/dobbinke/Desktop/Git/SSpred2/Programs/summer18redo/SampSize/Functions4sampsizeandtol.R")
source("~/Desktop/Git/SSpred2/Programs/summer18redo/SampSize/Functions4sampsizeandtol.R")

# file for testing functions in functions4sampsizeandtol.R

# test function just created
myparams <- get.params4tolerance(
	OptAcc = 0.8415,
	MCruns = 100,
	p1 = 0.5,
	NmbrInfFeatures = 9,
	CovMatStruct = "AR1.block",
	Block.size = 3,
	CorrParamVal = 0.7,
	Dimension = 500,
	trim.prop = 0.0,
	type.measure.for.cvglmnet = "class",
	number.of.computer.clusters = 11)
myparams;

# test function just created
mycovmat <- get.covmat(
	NmbrInfFeatures = myparams$NmbrInfFeatures,
	CovMatStruct = myparams$CovMatStruct,
	Block.size = myparams$Block.size,
	CorrParamVal = myparams$CorrParamVal,
	Dimension = myparams$Dimension) 
dim(mycovmat)
mycovmat[1:10,1:10]
covmat.Inverse = matrix.inverse(mycovmat)


# test function just created
mulength.TransformedSpace = get.mulength.TransformedSpace(
	p1 = myparams$p1,
	OptAcc = myparams$OptAcc)
mulength.TransformedSpace

my.meanvec.untransformed <- get.meanvec.UntransformedSpace(
	NmbrInfFeatures = myparams$NmbrInfFeatures,
	mulength.TransformedSpace = mulength.TransformedSpace,
	covmat.Inverse = covmat.Inverse) 
my.meanvec.untransformed
# Note that vector length is Dimension

# test function just created
install.packages("expm")
library(expm)
covmat.Root = sqrtm(mycovmat)
covmat.Inverse = matrix.inverse(mycovmat)
covmat.Inverse.Root <- sqrtm(covmat.Inverse)
covmat.InfFeatures <- mycovmat[1:myparams$NmbrInfFeatures,1:myparams$NmbrInfFeatures]
my.meanvec.transformed <- get.meanvec.TransformedSpace(
	my.meanvec.untransformed,
	covmat.Inverse.Root,
	myparams$NmbrInfFeatures)
my.meanvec.transformed

# test function just created
get.betainfty(my.meanvec.transformed,covmat.Inverse);

# # test function just created
# covmat.Root <- sqrtm(mycovmat)
# tolerance.type = "Logistic"
# current.N = 400;
# thistol <- EstimateTol(
# 	current.N = current.N,
# 	MCruns = myparams$MCruns,
#     covmat = mycovmat,
#     Dimension = myparams$Dimension,
#     covmat.Root = sqrtm(mycovmat),
#     type.measure = myparams$type.measure,
#     number.of.computer.clusters = myparams$number.of.computer.clusters,
#     p1 = myparams$p1,
# 	meanvec.UntransformedSpace = my.meanvec.untransformed,
# 	AccuracyORLogistic = "Logistic",
# 	NmbrInfFeatures = myparams$NmbrInfFeatures, 
#     OptAcc = myparams$OptAcc)
# thistol

# tolerance.type = "Accuracy"
# thistol <- EstimateTol(
# 	current.N = current.N,
# 	MCruns = myparams$MCruns,
#     covmat = mycovmat,
#     Dimension = myparams$Dimension,
#     covmat.Root = sqrtm(mycovmat),
#     type.measure = myparams$type.measure,
#     number.of.computer.clusters = myparams$number.of.computer.clusters,
#     p1 = myparams$p1,
# 	meanvec.UntransformedSpace = my.meanvec.untransformed,
# 	AccuracyORLogistic = "Logistic",
# 	NmbrInfFeatures = myparams$NmbrInfFeatures, 
#     OptAcc = myparams$OptAcc)
# thistol


# NOW TO TEST THE SAMPLE SIZE FUNCTION
tolerance.type <- 
AccuracyORLogistic <- "Logistic"
myparamsSampleSize <- 
  list(
  	TolTarget = 0.02,
  	LB.N = 1100,
  	UB.N = 2000,
  	meanvec.UntransformedSpace = my.meanvec.untransformed,
  	ToleranceType = tolerance.type,
  	MCruns = myparams$MCruns,
  	maxit = 20,
    covmat = mycovmat,
    Dimension = myparams$Dimension,
    covmat.Root = covmat.Root,
    type.measure = myparams$type.measure,
    number.of.computer.clusters = myparams$number.of.computer.clusters,
    p1 = myparams$p1,
    AccuracyORLogistic = AccuracyORLogistic,
    NmbrInfFeatures = myparams$NmbrInfFeatures,
    OptAcc = myparams$OptAcc
  	)
names(myparamsSampleSize)

oldwarn = getOption("warn")
options(warn = -1)
begin.time <- Sys.time();
SampSizeFun(
	TolTarget = myparamsSampleSize$TolTarget,
  	LB.N = myparamsSampleSize$LB.N,
  	UB.N = myparamsSampleSize$UB.N,
  	meanvec.UntransformedSpace = myparamsSampleSize$meanvec.UntransformedSpace,
  	MCruns = myparamsSampleSize$MCruns,
  	maxit = myparamsSampleSize$maxit,
    covmat = myparamsSampleSize$covmat ,
    Dimension = myparamsSampleSize$Dimension ,
    covmat.Root = myparamsSampleSize$covmat.Root ,
    type.measure = myparamsSampleSize$type.measure ,
    number.of.computer.clusters = myparamsSampleSize$number.of.computer.clusters ,
    p1 = myparamsSampleSize$p1 ,
    AccuracyORLogistic = myparamsSampleSize$AccuracyORLogistic ,
    NmbrInfFeatures = myparamsSampleSize$NmbrInfFeatures ,
    OptAcc = myparamsSampleSize$OptAcc 
  )
end.time <- Sys.time();
options(warn = oldwarn)
print(end.time - begin.time);

print(myresultsTable)

last.iter=myparamsSampleSize$maxit;
for (i in myparamsSampleSize$maxit:1)  {  # by reversing counting order, we avoid a need to break
  if (myresults$Table$tolest.path[i] == "NA") {
    last.iter <- (i-1);
  }
}


plot(1:last.iter,
    myresultsTable$tolest.path[1:last.iter],
    xlab="Step",
    ylab="Tol. Est.",
    type="lines")

for (i in 1:last.iter)  {
    text(i,myresultsTable$tolest.path[i],
        labels=as.character(myresultsTable$N.path[i]))
}


beep("fanfare")



