library(spBayes)

rf.n200.dat <- read.table("rf.n200.dat", header=TRUE)
rf.n200.dat <- as.data.frame(rf.n200.dat)


Y <- rf.n200.dat$Y
coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
w <- rf.n200.dat$w

##############################
##Simple spatial regression
##############################
m.1 <- spLM(Y~1, coords=coords,
             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
               "tau.sq.IG"=c(2, 1)),
             cov.model="exponential",
             n.samples=1000, verbose=TRUE, n.report=100)
             
# export Y. In this case X is just intercept, which is also the only case that we need to consider here
write.table(Y, "Ydat.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)