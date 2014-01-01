# # # # # # # # # # # #
# Jiehua Chen
# Oct. 15th, 2012
# jc3288@columbia.edu
# # # # # # # # # # # #

#read in data
setwd("../data")
MIRdata <- read.csv("AfSIS-Core MIR first derivative.csv")
MIRdata$PlotCultMgd <- ifelse(c(MIRdata$PlotCultMgd)==1, NA, MIRdata$PlotCultMgd)
MIRdata$PlotCultMgd <- MIRdata$PlotCultMgd -2
MIRdata$Depth.std <- c(MIRdata$Depth.std)

#study of Total.Carbon

bart.MIRdata <- cbind(MIRdata[, c("Total.Carbon")], MIRdata[,(129:1876)])
bart.MIRdata.narm <- na.omit(bart.MIRdata)

y <- bart.MIRdata.narm[,1]
x <- as.matrix(bart.MIRdata.narm[,-1])

# load revised BART C++ code
setwd("../src")
dyn.load("mbart.so")
source("bartfunc.R")
source("makeind.R")

#set up folder to save the MCMC trees. You need to make sure that the folder is empty before running the following code.

bart.est <- bart_saveresults(x, y, sigest = sd(y), ndpost=50, nskip=500, keepevery=10)


