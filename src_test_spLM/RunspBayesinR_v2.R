library(spBayes)

rmvn <- function(n, mu=0, V = matrix(1)){
		p <- length(mu)
				if(any(is.na(match(dim(V),p))))
						stop("Dimension problem!")
								D <- chol(V)
								t(matrix(rnorm(n*p), ncol=p)%*%D + rep(mu,rep(n,p)))
}
set.seed(1)
n <- 100
coords <- cbind(runif(n,0,1), runif(n,0,1))
X <- as.matrix(cbind(1, rnorm(n)))
B <- as.matrix(c(1,5))
p <- length(B)
sigma.sq <- 2
tau.sq <- 0.1
phi <- 3/0.5
D <- as.matrix(dist(coords))
R <- exp(-phi*D)
w <- rmvn(1, rep(0,n), sigma.sq*R)
y <- rnorm(n, X%*%B + w, sqrt(tau.sq))
n.samples <- 2000
starting <- list("phi"=3/0.5, "sigma.sq"=50, "tau.sq"=1)
tuning <- list("phi"=0.1, "sigma.sq"=0.1, "tau.sq"=0.1)
priors.1 <- list("beta.Norm"=list(rep(0,p), diag(1000,p)),
				"phi.Unif"=c(3/1, 3/0.1), "sigma.sq.IG"=c(2, 2),
				"tau.sq.IG"=c(2, 0.1))
priors.2 <- list("beta.Flat", "phi.Unif"=c(3/1, 3/0.1),
				"sigma.sq.IG"=c(2, 2), "tau.sq.IG"=c(2, 0.1))
cov.model <- "exponential"
n.report <- 500
verbose <- TRUE
m.1 <- spLM(y~X-1, coords=coords, starting=starting,
				tuning=tuning, priors=priors.1, cov.model=cov.model,
				n.samples=n.samples, verbose=verbose, n.report=n.report)
write.table(y, "Ydat.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
#
#rf.n200.dat <- read.table("rf.n200.dat", header=TRUE)
#rf.n200.dat <- as.data.frame(rf.n200.dat)
#
#
#Y <- rf.n200.dat$Y
#coords <- as.matrix(rf.n200.dat[,c("x.coords","y.coords")])
#w <- rf.n200.dat$w
#
###############################
###Simple spatial regression
###############################
#m.1 <- spLM(Y~1, coords=coords,
#             starting=list("phi"=0.6,"sigma.sq"=1, "tau.sq"=1),
#             sp.tuning=list("phi"=0.01, "sigma.sq"=0.05, "tau.sq"=0.05),
#             priors=list("phi.Unif"=c(0.3, 3), "sigma.sq.IG"=c(2, 1),
#               "tau.sq.IG"=c(2, 1)),
#             cov.model="exponential",
#             n.samples=1000, verbose=TRUE, n.report=100)
#             
## export Y. In this case X is just intercept, which is also the only case that we need to consider here
#write.table(Y, "Ydat.txt", col.names=FALSE, row.names=FALSE, quote=FALSE)
