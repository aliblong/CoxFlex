library(MASS)
library(PermAlgo)
library(survival)



N<-300
t<-12
incidence.rate<-0.67
beta2<-log(1.5)
beta3<-log(1.5)
Sigma <- matrix(c(1,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,0,0,0,
                  0.7,1,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,0,0,
                  0.6,0.7,1,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,0,
                  0.5,0.6,0.7,1,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0,
                  0.4,0.5,0.6,0.7,1,0.7,0.6,0.5,0.4,0.3,0.2,0.1,
                  0.3,0.4,0.5,0.6,0.7,1,0.7,0.6,0.5,0.4,0.3,0.2,
                  0.2,0.3,0.4,0.5,0.6,0.7,1,0.7,0.6,0.5,0.4,0.3,
                  0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,0.7,0.6,0.5,0.4,
                  0  ,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,0.7,0.6,0.5,
                  0  ,0  ,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,0.7,0.6,
                  0  ,0  ,0  ,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1,0.7,
                  0  ,0  ,0  ,0  ,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1),12,12, byrow=T)
X1<-as.data.frame(mvrnorm(N, c(rep(0,t)), Sigma))
x1.long <- reshape(X1,times = names(X1), varying = list(names(X1)),direction = "long")
x1.long<-x1.long[order(x1.long$id),]
x2<-rnorm(N, 0, 1)
x3<-rep(exp(rnorm(N,0,1)), each=t)
g1<-2*exp(-(c(1:12)-6)^{2}/8)
gx<-2+(x1.long$V1)^3/8
g1x1<-g1*gx

Xmat<-NULL
Xmat<-cbind(x1.long$V1,g1x1,x2,x3)

eventtimes <- ceiling(rexp(N, -log(1-incidence.rate)/t))

censortimes <- rep(t, times=N)

sampledat<- permalgorithm(N, t, Xmat, XmatNames=c("x1","g1gx","x2","x3"),
                          eventRandom = eventtimes, censorRandom=censortimes,
                          betas=c(0,1,beta2,beta3), groupByD=FALSE )
sampledata<-sampledat[,-7]

devtools::use_data(sampledata, overwrite = TRUE)
