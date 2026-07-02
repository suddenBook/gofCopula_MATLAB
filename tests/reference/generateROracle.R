#!/usr/bin/env Rscript
# Generate deterministic reference values from gofCopula 0.4-3 and copula.

args <- commandArgs(trailingOnly = TRUE)
output <- if (length(args)) args[[1]] else "tests/reference/rOracle.json"
library(copula)
library(gofCopula)
library(jsonlite)

U <- matrix(c(
  .19,.31,
  .46,.73,
  .81,.58), ncol = 2, byrow = TRUE)

models <- list(
  normal = normalCopula(.35),
  t = tCopula(.25, df = 5, df.fixed = TRUE),
  clayton = claytonCopula(1.4),
  gumbel = gumbelCopula(1.7),
  frank = frankCopula(2.2),
  joe = joeCopula(1.8),
  amh = amhCopula(.45),
  galambos = galambosCopula(1.1),
  huslerreiss = huslerReissCopula(1.2),
  tawn = tawnCopula(.7),
  tev = tevCopula(.3, df = 5, df.fixed = TRUE),
  fgm = fgmCopula(-.6),
  plackett = plackettCopula(2.5))

theta <- c(normal=.35,t=.25,clayton=1.4,gumbel=1.7,frank=2.2,joe=1.8,
  amh=.45,galambos=1.1,huslerreiss=1.2,tawn=.7,tev=.3,fgm=-.6,plackett=2.5)
df <- c(normal=4,t=5,clayton=4,gumbel=4,frank=4,joe=4,amh=4,
  galambos=4,huslerreiss=4,tawn=4,tev=5,fgm=4,plackett=4)

familyValues <- lapply(names(models), function(name) {
  model <- models[[name]]
  list(theta=unname(theta[[name]]), df=unname(df[[name]]),
       cdf=as.numeric(pCopula(U,model)), pdf=as.numeric(dCopula(U,model)))
})
names(familyValues) <- names(models)

standard <- c("normal","t","clayton","gumbel","frank","joe","amh")
rosenblatt <- lapply(standard, function(name) unname(cCopula(U,models[[name]])))
names(rosenblatt) <- standard
rosenblatt$galambos <- unname(gofCopula:::.rosenblatt.galambos(U,theta[["galambos"]]))
rosenblatt$fgm <- unname(gofCopula:::.rosenblatt.fgm(U,theta[["fgm"]]))
rosenblatt$plackett <- unname(gofCopula:::.rosenblatt.plackett(U,theta[["plackett"]]))

S <- matrix(c(
  .08,.16,
  .18,.42,
  .31,.27,
  .44,.69,
  .57,.51,
  .68,.83,
  .79,.63,
  .92,.88), ncol=2, byrow=TRUE)
clayton <- claytonCopula(1.4)
Z <- cCopula(S,clayton)
AZ <- gofCopula:::.ArchmRtrans(clayton,S)
set.seed(92817)
referenceSample <- rCopula(500,clayton)
referenceValues <- copula:::F.n(referenceSample,referenceSample)

statistics <- list(
  cvm=unname(gofCopula:::.Tstats(S,"Sn",clayton)),
  ks=unname(gofCopula:::.Tstats(S,"KS",clayton)),
  kendallCvm=unname(gofCopula:::.Tstats(S,"SnK",clayton,cop.compare=referenceValues)),
  kendallKs=unname(gofCopula:::.Tstats(S,"TnK",clayton,cop.compare=referenceValues)),
  rosenblattSnB=unname(gofCopula:::.Tstats(Z,"SnB",clayton)),
  rosenblattSnC=unname(gofCopula:::.Tstats(Z,"SnC",clayton)),
  rosenblattGamma=unname(gofCopula:::.Tstats(Z,"AnGamma",clayton)),
  rosenblattChisq=unname(gofCopula:::.Tstats(Z,"AnChisq",clayton)),
  archmSnB=unname(gofCopula:::.Tstats(AZ,"SnB",clayton)),
  archmSnC=unname(gofCopula:::.Tstats(AZ,"SnC",clayton)),
  archmGamma=unname(gofCopula:::.Tstats(AZ,"AnGamma",clayton)),
  archmChisq=unname(gofCopula:::.Tstats(AZ,"AnChisq",clayton)))

# A moderate-dependence sample avoids known convergence failures in R's
# leave-one-out optimizer and isolates the PIOS formula comparison.
set.seed(321)
piosSample <- rCopula(12,clayton)
statistics$piosRn <- unname(gofCopula:::.Tstats(piosSample,"Rn",clayton))
statistics$piosTn <- unname(gofCopula:::.Tstats(piosSample,"Tn",clayton,
  add.parameters=list(nrow(piosSample),1,"mpl")))
statistics$white <- unname(VineCopula::BiCopGofTest(
  S[,1],S[,2],family=3,par=1.4,method="white",B=0)$statistic)
statistics$whiteT <- unname(VineCopula::BiCopGofTest(
  S[,1],S[,2],family=2,par=.25,par2=5,method="white",B=0)$statistic)

estimationModels <- list(
  normal=normalCopula(dim=2),
  t=tCopula(dim=2,df=5,df.fixed=TRUE),
  clayton=claytonCopula(dim=2),
  gumbel=gumbelCopula(dim=2),
  frank=frankCopula(dim=2),
  joe=joeCopula(dim=2))
estimates <- lapply(estimationModels, function(model) unname(
  fitCopula(model,S,method="mpl",estimate.variance=FALSE)@estimate))

set.seed(4815)
kernelSample <- rCopula(100,clayton)
kernelGrid <- SparseGrid::createIntegrationGrid("GQU",dimension=2,k=6)
kernelBandwidth <- as.vector(diag(2.6073*nrow(S)^(-1/6)*chol(cov(S)))*.5)
statistics$kernel <- sum(vapply(seq_len(nrow(kernelGrid$nodes)), function(i)
  kernelGrid$weights[i]*gofCopula:::.integrand(
    kernelGrid$nodes[i,],S,kernelSample,kernelBandwidth),numeric(1)))

set.seed(73421)
bootstrapCount <- 7
bootstrapFit <- fitCopula(claytonCopula(dim=2),piosSample,method="mpl",
  estimate.variance=FALSE)
bootstrapModel <- claytonCopula(bootstrapFit@estimate)
bootstrapObserved <- gofCopula:::.Tstats(piosSample,"Sn",bootstrapModel)
bootstrapSamples <- array(0,dim=c(nrow(piosSample),2,bootstrapCount))
bootstrapStatistics <- numeric(bootstrapCount)
for (b in seq_len(bootstrapCount)) {
  bootstrapSamples[,,b] <- rCopula(nrow(piosSample),bootstrapModel)
  replicateFit <- fitCopula(claytonCopula(dim=2),bootstrapSamples[,,b],
    method="mpl",estimate.variance=FALSE)
  bootstrapStatistics[b] <- gofCopula:::.Tstats(bootstrapSamples[,,b],"Sn",
    claytonCopula(replicateFit@estimate))
}
bootstrap <- list(samples=unname(bootstrapSamples),
  observed=unname(bootstrapObserved),statistics=unname(bootstrapStatistics),
  pValue=mean(abs(bootstrapStatistics)>=abs(bootstrapObserved)))

# --- Trivariate cases -------------------------------------------------------
U3 <- matrix(c(
  .22,.41,.36,
  .48,.67,.55,
  .71,.29,.62,
  .87,.79,.91), ncol=3, byrow=TRUE)
d3models <- list(
  normal=normalCopula(.45, dim=3, dispstr="ex"),
  clayton=claytonCopula(1.6, dim=3),
  gumbel=gumbelCopula(1.9, dim=3),
  frank=frankCopula(3.1, dim=3),
  joe=joeCopula(2.1, dim=3))
d3theta <- c(normal=.45, clayton=1.6, gumbel=1.9, frank=3.1, joe=2.1)
set.seed(55611)
S3 <- rCopula(10, claytonCopula(1.6, dim=3))
set.seed(55612)
E3 <- rCopula(12, claytonCopula(1.6, dim=3))
d3 <- lapply(names(d3models), function(name) {
  model <- d3models[[name]]
  Z3 <- cCopula(U3, model)
  ZS <- cCopula(S3, model)
  list(theta=unname(d3theta[[name]]),
       cdf=as.numeric(pCopula(U3, model)),
       pdf=as.numeric(dCopula(U3, model)),
       rosenblatt=unname(Z3),
       cvm=unname(gofCopula:::.Tstats(S3,"Sn",model)),
       ks=unname(gofCopula:::.Tstats(S3,"KS",model)),
       rosenblattSnB=unname(gofCopula:::.Tstats(ZS,"SnB",model)),
       rosenblattSnC=unname(gofCopula:::.Tstats(ZS,"SnC",model)),
       rosenblattGamma=unname(gofCopula:::.Tstats(ZS,"AnGamma",model)),
       rosenblattChisq=unname(gofCopula:::.Tstats(ZS,"AnChisq",model)),
       estimate=unname(fitCopula(
         switch(name,
           normal=normalCopula(dim=3, dispstr="ex"),
           clayton=claytonCopula(dim=3),
           gumbel=gumbelCopula(dim=3),
           frank=frankCopula(dim=3),
           joe=joeCopula(dim=3)),
         E3, method="mpl", estimate.variance=FALSE)@estimate))
})
names(d3) <- names(d3models)

# --- Rank-margins bootstrap pipeline (rCompatible parity anchor) -------------
set.seed(88231)
rawX <- cbind(rnorm(12, 3, 2), rexp(12, .8))
rankedX <- apply(rawX, 2, function(col) ecdf(col)(col) * 12 / 13)
ranksFit <- fitCopula(claytonCopula(dim=2), rankedX, method="mpl",
  estimate.variance=FALSE)
ranksModel <- claytonCopula(ranksFit@estimate)
ranksObserved <- gofCopula:::.Tstats(rankedX, "Sn", ranksModel)
set.seed(88232)
ranksCount <- 7
ranksSamples <- array(0, dim=c(12, 2, ranksCount))
ranksStatistics <- numeric(ranksCount)
for (b in seq_len(ranksCount)) {
  ranksSamples[,,b] <- rCopula(12, ranksModel)
  replicateFit <- fitCopula(claytonCopula(dim=2), ranksSamples[,,b],
    method="mpl", estimate.variance=FALSE)
  # R computes replicate statistics on the RAW draws (no re-ranking).
  ranksStatistics[b] <- gofCopula:::.Tstats(ranksSamples[,,b], "Sn",
    claytonCopula(replicateFit@estimate))
}
ranksPipeline <- list(rawData=unname(rawX), ranked=unname(rankedX),
  theta=unname(ranksFit@estimate), observed=unname(ranksObserved),
  samples=unname(ranksSamples), statistics=unname(ranksStatistics),
  pValue=mean(abs(ranksStatistics) >= abs(ranksObserved)))

# --- White test, Student copula, moderate sample size ------------------------
set.seed(99417)
whiteSample <- rCopula(100, tCopula(.4, df=6, df.fixed=TRUE))
whiteT100 <- list(sample=unname(whiteSample), par=.4, df=6,
  statistic=unname(VineCopula::BiCopGofTest(whiteSample[,1], whiteSample[,2],
    family=2, par=.4, par2=6, method="white", B=0)$statistic))

# --- Parametric margin fits ---------------------------------------------------
set.seed(77123)
marginData <- list(
  norm=rnorm(30, 5, 2),
  exp=rexp(30, .7),
  gamma=rgamma(30, 2.5, 1.3),
  lnorm=rlnorm(30, .8, .5),
  weibull=rweibull(30, 1.8, 2.2),
  beta=rbeta(30, 2, 3),
  cauchy=rcauchy(30, 1, .8),
  f=rf(30, 8, 12))
margins <- lapply(names(marginData), function(name) {
  fit <- gofCopula:::.one.mar(marginData[[name]], name)
  list(data=marginData[[name]], parameters=as.numeric(fit[[1]]),
       u=as.numeric(fit[[2]]))
})
names(margins) <- names(marginData)

# --- Flip conventions (Frank admits every rotation) ---------------------------
flip <- lapply(c(90, 180, 270), function(angle) {
  xf <- gofCopula:::.rotateCopula(S, angle)
  fit <- fitCopula(frankCopula(dim=2), xf, method="mpl",
    estimate.variance=FALSE)
  list(angle=angle, theta=unname(fit@estimate),
       cvm=unname(gofCopula:::.Tstats(xf, "Sn", frankCopula(fit@estimate))))
})

# --- t-EV joint (theta, df) estimation ----------------------------------------
set.seed(66101)
tevSample <- rCopula(100, tevCopula(.5, df=4, df.fixed=TRUE))
tevFit <- fitCopula(tevCopula(df.fixed=FALSE), tevSample, method="mpl",
  estimate.variance=FALSE)
tev <- list(sample=unname(tevSample),
  theta=unname(tevFit@estimate[1]), df=unname(tevFit@estimate[2]))

oracle <- list(
  metadata=list(gofCopula=as.character(packageVersion("gofCopula")),
    copula=as.character(packageVersion("copula")),
    VineCopula=as.character(packageVersion("VineCopula")),
    R=as.character(getRversion())),
  evaluationPoints=unname(U), families=familyValues,
  rosenblatt=rosenblatt, statisticSample=unname(S),
  piosSample=unname(piosSample), kernelSample=unname(kernelSample),
  kendallReference=as.numeric(referenceValues),
  statistics=statistics, estimates=estimates, bootstrap=bootstrap,
  evaluationPoints3=unname(U3), statisticSample3=unname(S3),
  estimationSample3=unname(E3), d3=d3,
  ranksPipeline=ranksPipeline, whiteT100=whiteT100,
  margins=margins, flip=flip, tev=tev)

write_json(oracle,output,auto_unbox=TRUE,digits=17,pretty=TRUE)
