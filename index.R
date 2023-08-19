#| include: false
library(np)
library(crs)
library(car)
library(KernSmooth)
library(latex2exp)
library(MASS)
library(lmtest)
options(crs.messages=FALSE,np.messages=FALSE,np.tree=TRUE)

#| echo: true
#| eval: false
## ## Generate some data: sex (unordered categorical), income (ordered categorical),
## ## and height (numeric)
## n <- 100
## sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
## income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
## height <- rnorm(n,mean=150,sd=20)
## ## Note - by default these variables may not be the data types we desire
## class(sex);class(income);class(height)
## ## income is already numeric(), but sex and height are character()
## ## sex is categorical and unordered, so cast as type factor()
## sex <- factor(sex)
## ## income is categorical and ordered, but we need to ensure intended order (it
## ## will assume alphabetical ordering otherwise). Suppose you ignore it - let's
## ## see what happens when we just cast as type ordered() using defaults
## income <- ordered(income)
## ## The levels are in alphabetical order, which we don't want
## levels(income)
## ## We shall reorder the ordered factor levels as intended using levels=...
## income <- ordered(income,levels=c("Low","Middle","High"))
## levels(income)
## ## Check data types again
## class(sex);class(income);class(height)
## ## Note that with integers the default ordered() works fine
## x <- sample(c(2,5,4,3,1),n,replace=TRUE)
## x <- ordered(x)
## levels(x)


#| eval: true
#| echo: true
## Let's simulate a random numeric sample from the normal distribution
set.seed(42)
n <- 1000
x <- rnorm(n)
## Let's sort the data so we can graph x versus dnorm(x,...) using lines (type="l")
x <- sort(x)
## Since we simulated the data, let's plot the true, known, parametric density
## (we can't do this with actual data because the density of such data is, in
## general, unknown)
plot(x,dnorm(x,mean=mean(x),sd=sd(x)),type="l",ylab="Parametric Density Estimate",xlab="X")



pander::pander(shapiro.test(x))



data(faithful)
with(faithful,pander::pander(shapiro.test(eruptions)))



plot(density(faithful$eruptions),main="")
rug(faithful$eruptions)



plot(density(faithful$eruptions),main="")
with(faithful,lines(density(eruptions)$x,dnorm(density(eruptions)$x,mean=mean(eruptions),sd=sd(eruptions)),col=2,lty=2))
rug(faithful$eruptions)
legend("topleft",c("Nonparametric","Parametric"),lty=c(1,2),col=c(1,2),bty="n")



hist(faithful$eruptions,prob=TRUE,main="",xlab="Eruptions",breaks=20,xlim=c(1.25,5.5))
with(faithful,lines(density(eruptions)$x,fitted(npudens(tdat=eruptions,edat=density(eruptions)$x))))
rug(faithful$eruptions)


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions,data=faithful)
plot(fhat,neval=250,plot.errors.method="bootstrap")


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions+waiting,data=faithful)
plot(fhat,theta=330,xtrim=-0.05,neval=75,view="fixed",main="")


#| echo: true
n <- 250
set.seed(42)
sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
sex <- factor(sex)
phat <- npudens(~sex)
plot(phat,plot.errors.method="bootstrap")


#| echo: true
n <- 250
set.seed(42)
income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
income <- ordered(income,levels=c("Low","Middle","High"))
phat <- npudens(~income)
plot(phat,plot.errors.method="bootstrap")


#| label: tbl-wage1mixedtable
#| tbl-cap: Counts of number of dependants present in 526 households by cell
library(np)
library(plot3D)
data(wage1)
knitr::kable(with(wage1,t(data.frame(numdep=sort(unique(numdep)),counts=as.numeric(table(numdep))))),
             booktabs=TRUE,
             linesep="")


#| echo: true
#| eval: false
## library(np)
## library(plot3D)
## data(wage1)
## numdep.seq <- with(wage1,sort(unique(numdep)))
## lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
## wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
## bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
## fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## ## Hack since scatter3D converts ordered 0-6 to numeric 1-7
## scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
##           ylab="Log-Wage",
##           xlab="Number of Dependants",
##           zlab="Joint Density",
##           ticktype="detailed",
##           angle=15,
##           box=TRUE,
##           type="h",
##           grid=TRUE,
##           col="blue",
##           colkey=FALSE)


#| label: fig-wage1mixeddensity
#| fig-cap: Mixed-data bivariate kernel density estimate for the joint PDF of lwage (numeric) and numdeps (ordered)
data(wage1)
bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
numdep.seq <- with(wage1,sort(unique(numdep)))
lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## Hack since scatter3D converts ordered 0-6 to numeric 1-7
scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
          ylab="Log-Wage",
          xlab="Number of Dependants",
          zlab="Joint Density",
          ticktype="detailed",
          angle=15,
          box=TRUE,
          type="h",
          grid=TRUE,
          col="blue",
          colkey=FALSE)


rm(list=ls())
library(np)
library(car)
Prestige <- Prestige[order(Prestige$income),]
attach(Prestige)
## Plug-in
library(KernSmooth)
h <- dpill(income, prestige)
model.plugin <- locpoly(income, prestige, bandwidth = h)
## Under
bw <- npregbw(prestige~income,bws=c(sd(income)*length(income)^{-1/5}/10),regtype="ll",bandwidth.compute=FALSE)
fit.under <- fitted(npreg(bws=bw))
## Over
bw <- npregbw(prestige~income,bws=c(1000*sd(income)*length(income)^{-1/5}),regtype="ll",bandwidth.compute=FALSE)
fit.over <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.ls",regtype="ll")
fit.ls <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.aic",regtype="ll")
fit.aic <- fitted(npreg(bws=bw))
par(mfrow=c(2,2))
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Undersmoothed")
lines(income,fit.under)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Oversmoothed")
lines(income,fit.over)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Plug-In")
lines(model.plugin$x,model.plugin$y)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="AICc & CV-LS")
lines(income,fit.ls,lty=1)
lines(income,fit.aic,lty=2)


#| echo: true
library(np)
set.seed(42)
n <- 1000
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + rnorm(n,sd=0.25*sd(dgp))
ghat <- npreg(y~x)
plot(x,y,cex=0.5,col="grey",xlab="X",ylab="Y")
lines(x,dgp)
lines(x,fitted(ghat),col=2)
abline(ghat.ols <- lm(y~x),col=3)
legend("top",c("DGP","Kernel","OLS"),col=1:3,lty=1,bty="n")



pander::pander(summary(ghat.ols))


#| echo: true
plot(ghat,gradients=TRUE,neval=250)
lines(x,-2*pi*sin(2*pi*x),col=2)
abline(h=coef(ghat.ols)[2],col=3)
legend("topleft",c("Kernel ME","DGP ME","Linear ME"),col=1:3,lty=1,bty="n")


#| echo: true
library(np)
data(wage1)
## For details see ?wage1
ghat <- npreg(lwage ~ female + married + educ + exper + tenure, data=wage1, regtype="ll")
summary(ghat)


#| echo: true
## We run out of graph axis dimensions with > 2 predictors, so it is common to
## construct partial plots that plot the fitted model versus each predictor
## separately holding the off-axis predictors at, say, their median value (you
## can change this - see ?npplot and the argument xq)
par(mfrow=c(2,3))
plot(ghat,plot.errors.method="bootstrap")


#| echo: true
par(mfrow=c(2,3))
plot(ghat,gradients=TRUE,plot.errors.method="bootstrap")


#| echo: true
ghat.ols <- lm(lwage ~ female + married + educ + exper + tenure, data=wage1)
summary(ghat.ols)


#| echo: true
npsigtest(ghat)


#| echo: true
attach(wage1)
df <- data.frame(female = factor("Male", levels=levels(female)),
                 married = factor("Notmarried", levels=levels(married)),
                 educ = median(educ),
                 tenure = median(tenure),
                 exper = median(exper))
head(df)
predict(ghat, newdata=df)
## Or you could use ghat <- npreg(...,newdata=df) and fitted(ghat) 
## ghat <- npreg(lwage ~ female + married + educ + exper + tenure, 
##               data=wage1,
##               regtype="ll", 
##               newdata=df)
## fitted(ghat)
## Note - If `df` had multiple rows we would get a prediction for each
## row (i.e., a vector of predictions)

#| echo: true
#| eval: false
## ## Example 1 - create a simple time series that begins on the 2nd Quarter of
## ## 1959
## x <- ts(1:10, frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)
## ## Example 2 - create simulated data from an ARIMA model that begins on the 2nd
## ## Quarter of 1959
## x <- ts(arima.sim(n = 10, list(ar = c(0.8897, -0.4858),
##                                ma = c(-0.2279, 0.2488)),
##                   sd = sqrt(0.1796)),
##         frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)


#| echo: true
#| eval: false
## ## ldeaths is monthly deaths from bronchitis, emphysema and asthma in the UK,
## ## 1974–1979, both sexes
## ldeaths
## ## lag(x, k=1, ...), note the sign of k: a series lagged by a positive k, i.e.,
## ## y_{t+k}
## lag(ldeaths, 0) # returns original series y_t
## lag(ldeaths, -1) # returns series lagged once y_{t-1}
## lag(ldeaths, -2) # returns series lagged twice y_{t-2}
## ## In matrix form compare to cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))
## cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))


#| echo: true
#| eval: false
## diff(ldeaths)
## ## Manually compute first difference
## ldeaths[2:length(ldeaths)]-ldeaths[1:(length(ldeaths)-1)]


#| echo: true
#| eval: false
## library(np)
## ## Simulate a simple AR(1) process with lag coefficient 0.9
## set.seed(42)
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Using lag() in lm() calls will fail (same result regardless of lag, R^2=1,
## ## lag() is totally ignored so you are regressing y_t on y_t though naturally
## ## you may not suspect that is what would happen in this case)
## ghat.ols <- lm(y~lag(y,-1))
## summary(ghat.ols)
## ## But if you manually regress numeric y_t on y_{t-1} you get the expected
## ## result (try also ar.ols(y, order.max = 1) or arima(y, order = c(1,0,0))).
## ## Here we regress (y_t,...,y_2) on (y_{t-1},...,y_1) using numeric vectors
## ghat.ols <- lm(y[2:length(y)]~y[1:(length(y)-1)])
## summary(ghat.ols)
## ## npreg() does support ts() objects. In the example below we use the npreg()
## ## function to conduct simple OLS by noting that local linear regression with a
## ## large bandwidth is linear least squares, so the gradient will be a vector of
## ## constants equal to the least squares coefficient for lag(y,-1) - think of the
## ## lag() function and notation y_{t-1} or "why tee minus one" hence the -1 here
## ghat <- npreg(y~lag(y,-1),regtype="ll",bws=10^5,bandwidth.compute=FALSE,gradients=TRUE)
## tail(gradients(ghat))


#| echo: true
#| eval: false
## library(np)
## ## Simulate data from a stationary univariate AR(1) process
## set.seed(42)
## n <- 100
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Conduct local constant estimation of y on four lags
## ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## ## Re-fit the model using only 1 lag
## ghat.1.lag <- npreg(y~lag(y,-1))
## ## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## ## Plot the data and nonparametric fit
## plot(y,lwd=2)
## lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
## points(ts(fitted(ghat.1.lag),start=2),col=3)
## legend("topright",c("Series",
##                     "npreg() including 4 lags",
##                     "npreg() including 1 lag"),
##                     lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")


#| echo: false
#| label: fig-arimatsnpreg
#| fig-cap: Local constant estimation of a stationary univariate AR(1) time series, $Y_t=0.9Y_{t-1}+\epsilon_t$
library(np)
## Simulate data from a stationary univariate AR(1) process
set.seed(42)
n <- 100
y <- arima.sim(n = n, list(ar = c(0.9)))
## Conduct local constant estimation of y on four lags
ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## Re-fit the model using only 1 lag
ghat.1.lag <- npreg(y~lag(y,-1))
## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## Plot the data and nonparametric fit
plot(y,lwd=2)
lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
points(ts(fitted(ghat.1.lag),start=2),col=3)
legend("topright",c("Series",
                    "npreg() including 4 lags",
                    "npreg() including 1 lag"),
                    lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")



## Examine the model summary
summary(ghat.4.lags)


#| echo: true
## Conduct a nonparametric significance test
npsigtest(ghat.4.lags)


#| echo: true
## npglpreg() does not support time series objects, so we manually
## create the lagged series
library(crs)
options(crs.messages=FALSE)
y.lag.1 <- c(rep(NA,1),y[1:(n-1)])
y.lag.2 <- c(rep(NA,2),y[1:(n-2)])
y.lag.3 <- c(rep(NA,3),y[1:(n-3)])
y.lag.4 <- c(rep(NA,4),y[1:(n-4)])
model.glp <- npglpreg(y~y.lag.1+y.lag.2+y.lag.3+y.lag.4,nmulti=5,degree.max=5)


#| echo: true
summary(model.glp)


#| echo: false
model.lm <- lm(y~y.lag.1,subset=5:n)


#| echo: true
summary(model.lm)


#| echo: true
## Examine the first few fitted values from the nonparametric and
## linear AR(1) models
cbind(fitted(model.glp),fitted(model.lm))[1:10,]


#| echo: true
#| eval: false
## ## Fan & Yao's (1998) estimator, trim negative values of they occur
## ## Chen, Cheng & Peng's (2009) estimator is always positive but my
## ## Monte Carlo simulations show it is less efficient for this DGP
## library(np)
## set.seed(42)
## n <- 1000
## x <- sort(runif(n))
## sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
## dgp <- sin(2*pi*x)
## y <- dgp + rnorm(n,sd=sqrt(sigma.var))
## model <- npreg(y~x,regtype="ll")
## r <- residuals(model)**2
## ## Fan and Yao's (1998) estimator with trimming if needed
## var.fy <- fitted(npreg(r~x,regtype="ll"))
## var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## ## Chen, Cheng and Peng's (2009) estimator
## log.r <- log(r+1/n) ## Avoids log(0)
## V.hat <- fitted(npreg(log.r~x,regtype="ll"))
## d.hat <- 1/mean(r*exp(-V.hat))
## var.ccp <- exp(V.hat)/d.hat
## par(mfrow=c(1,2))
## plot(x,y,cex=.25,col="grey")
## lines(x,dgp,col=1,lty=1,lwd=2)
## lines(x,fitted(model),col=2,lty=2,lwd=2)
## legend("topleft",c("g(x)=sin(2 pi x)",
##                    "Nonparametric Estimate"),
##        col=1:2,
##        lty=1:2,
##        lwd=c(2,2),
##        bty="n")
## ylim=c(min(r),quantile(r,0.95))
## plot(x,r,ylim=ylim,cex=.25,col="grey")
## lines(x,sigma.var,col=1,lty=1,lwd=2)
## lines(x,var.fy,col=2,lty=2,lwd=2)
## lines(x,var.ccp,col=3,lty=4,lwd=2)
## legend("topleft",c("volatility=(2 pi x)^2/10",
##                    "Fan and Yao (1998)",
##                    "Chen et al. (2009)"),
##        col=1:3,
##        lty=c(1,2,4),
##        lwd=c(2,2,2),
##        bty="n")


#| label: fig-fanyaosim
#| fig-cap: Fan and Yao's (1998) and Chen et al.'s (2009) conditional variance estimators
#| echo: false
par(mfrow=c(1,2))
## Fan & Yao's (1998) estimator, trim negative values of they occur
## Chen, Cheng & Peng's (2009) estimator is always positive but my
## Monte Carlo simulations show it is less efficient for this DGP
set.seed(42)
n <- 1000
x <- sort(runif(n))
sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
dgp <- sin(2*pi*x)
y <- dgp + rnorm(n,sd=sqrt(sigma.var))
model <- npreg(y~x,regtype="ll")
r <- residuals(model)**2
## Fan and Yao's (1998) estimator with trimming if needed
var.fy <- fitted(npreg(r~x,regtype="ll"))
var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## Chen, Cheng and Peng's (2009) estimator
log.r <- log(r+1/n) ## Avoids log(0)
V.hat <- fitted(npreg(log.r~x,regtype="ll"))
d.hat <- 1/mean(r*exp(-V.hat))
var.ccp <- exp(V.hat)/d.hat
plot(x,y,cex=.25,col="grey")
lines(x,dgp,col=1,lty=1,lwd=2)
lines(x,fitted(model),col=2,lty=2,lwd=2)
legend("topleft",c("g(x)=sin(2 pi x)", 
                   "Nonparametric Estimate"),
       col=1:2,
       lty=1:2,
       lwd=c(2,2),
       bty="n")
ylim=c(min(r),quantile(r,0.95))
plot(x,r,ylim=ylim,cex=.25,col="grey")
lines(x,sigma.var,col=1,lty=1,lwd=2)
lines(x,var.fy,col=2,lty=2,lwd=2)
lines(x,var.ccp,col=3,lty=4,lwd=2)
legend("topleft",c("volatility=(2 pi x)^2/10",
                   "Fan and Yao (1998)",
                   "Chen et al. (2009)"),
       col=1:3,
       lty=c(1,2,4),
       lwd=c(2,2,2),
       bty="n")


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## ## Compute the pooled first-step estimator
## model.pooled <- npreg(y~x,ckerorder=4)
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~x,ckerorder=2)


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## model.fe <- npreg(y~x+factor(id))


#| echo: true
#| eval: false
## library(Ecdat)
## data(Airline)
## ?Airline
## library(np)
## options(np.messages=FALSE)
## attach(Airline)
## 
## ## Dependent variable is (log) cost, predictors (log) output, fuel price, load factor, year
## lcost <- as.numeric(log(cost))
## loutput <- as.numeric(log(output))
## lpf <- as.numeric(log(pf))
## lf <- as.numeric(lf)
## year <- ordered(year)
## airline <- factor(airline)
## 
## ## Airline specific fixed effects
## 
## model.fe <- npreg(lcost~loutput+lpf+lf+year+airline,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine",
##                   gradients=TRUE)
## 
## summary(model.fe)
## summary(model.fe)
## ## Plot partial means
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)
## 
## ## Random effects
## 
## ## Compute the pooled first-step estimator
## 
## id <- airline
## n <- length(unique(airline))
## T <- length(unique(year))
## model.pooled <- npreg(lcost~loutput+lpf+lf+year,
##                       ckerorder=4,
##                       regtype="ll",
##                       bwmethod="cv.aic",
##                       ukertype="liracine",
##                       okertype="liracine")
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~loutput+lpf+lf+year,
##                   ckerorder=2,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine")
## 
## summary(model.re)
## ## Plot partial means
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)

#| echo: true
#| eval: true
library(np)
data(wage1)
model.lm <- lm(lwage~female+married+educ+tenure+exper+I(exper^2),data=wage1)
knitr::kable(coef(model.lm),col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.pl <- npplreg(lwage~female+married+educ+tenure|exper,
                    data=wage1)
knitr::kable(coef(model.pl),digits=7,col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.scoef <- npscoef(lwage~educ+tenure+exper+expersq|female+married,
                       data=wage1,
                       betas=TRUE)
summary(model.scoef)


#| echo: true
#| eval: true
## The smooth coefficients are vectors, one for each predictor in X (they "vary"
## with Z), so we compute the means of the columns of these coefficients, one for
## each X predictor
knitr::kable(colMeans(coef(model.scoef)),digits=7,col.names="Avg. Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.index <- npindex(lwage~educ+
                       tenure+
                       exper+
                       expersq+
                       female+
                       married,
                       method="ichimura",
                       data=wage1)
summary(model.index)


#| echo: true
#| eval: true
library(MASS)
data("birthwt")
model.logit <- glm(low~factor(smoke)+
                   factor(race)+
                   factor(ht)+
                   factor(ui)+
                   ordered(ftv)+
                   age+
                   lwt,
                   family=binomial(link=logit),data=birthwt)
cm.logit <- with(birthwt,table(low, ifelse(fitted(model.logit)>0.5, 1, 0)))
ccr.logit <- sum(diag(cm.logit))/sum(cm.logit)
knitr::kable(cm.logit)


#| echo: true
#| eval: true
library(MASS)
data(birthwt)
library(np)
model.index <- npindex(low~factor(smoke)+
                       factor(race)+
                       factor(ht)+
                       factor(ui)+
                       ordered(ftv)+
                       age+
                       lwt,
                       method="kleinspady",
                       data=birthwt)
cm.index <- with(birthwt,table(low, ifelse(fitted(model.index)>0.5, 1, 0)))
ccr.index <- sum(diag(cm.index))/sum(cm.index)
knitr::kable(cm.index)


require(np)
require(MASS)
## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(np)
options(np.messages=FALSE)
set.seed(42)
## This function takes a confusion matrix and formats it correctly if
## it is unbalanced and returns the CCR as well.
CM <- function(cm) {
  factor.values.eval <- colnames(cm)
  CM <- matrix(0,nrow(cm),nrow(cm))
  rownames(CM) <- rownames(cm)
  colnames(CM) <- rownames(cm)
  for(i in 1:ncol(cm)) CM[,(1:nrow(cm))[rownames(cm)==factor.values.eval[i]]] <- cm[,i]
  return(list(CM=CM,CCR=sum(diag(CM))/sum(CM)))
}
## Load the birthwt data
library(MASS)
data(birthwt)
## Create a data frame that has up to 4th order polynomials. We will
## use the BIC-optimal parametric model that will allow for
## interactions and some nonlinearity.
bwt <- with(birthwt,data.frame(low=factor(low),
                               race=factor(race),
                               smoke=factor(smoke),
                               ht=factor(ht),
                               ui=factor(ui),
                               ftv=ordered(ftv),
                               age=age,
                               agesq=age**2,
                               agecu=age**3,
                               agequ=age**4,
                               lwt=lwt,
                               lwtsq=lwt**2,
                               lwtcu=lwt**3,
                               lwtqu=lwt**4))
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), number of multistarts for bandwidth selection (nmulti)
## and number of train/eval splits (M)
n.eval <- 10
n <- nrow(bwt)
n.train <- n-n.eval
M <- 1000
## Create storage vectors
ccr.linear <- numeric(M)
ccr.linint <- numeric(M)
ccr.BIC <- numeric(M)
ccr.kernel <- numeric(M)
ccr.index <- numeric(M)
## Copy the full sample into the object train
train <- bwt
## Fit the parametric Logit model for the full sample.
##
## Linear 
logit.linear <- glm(low~
                    smoke+
                    race+
                    ht+
                    ui+
                    age+
                    lwt+
                    ftv,
                    family=binomial(link=logit),
                    data=train)
## Linear Logit model with interactions
logit.linint <- glm(low~
                    (smoke+
                     race+
                     ht+
                     ui+
                     age+
                     lwt+
                     ftv)^2,
                    family=binomial(link=logit),
                    data=train)
## BIC-optimal Logit model
logit.BIC <- glm(low ~ .,
                 family = binomial(link=logit),
                 data = train)
logit.BIC <- stepAIC(logit.BIC, ~ .^3,
                     trace=FALSE,
                     k=log(nrow(birthwt)))
## Klein spady estimator
model.index <- npindex(low~smoke+
                       race+
                       ht+
                       ui+
                       ftv+
                       age+
                       lwt,
                       method="kleinspady",
                       data=train)
## Get the bandwidths for the nonparametric model for the full sample.
bw <- npcdensbw(low~
                smoke+
                race+
                ht+
                ui+
                age+
                lwt+
                ftv,
                data=train)
## Apparent (in-sample) performance
ccr.app.linear <- with(train,CM(table(low,ifelse(predict(update(logit.linear),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.linint <- with(train,CM(table(low,ifelse(predict(update(logit.linint),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.BIC <- with(train,CM(table(low,ifelse(predict(update(logit.BIC),
                                              type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.index <- with(train,CM(table(low,ifelse(fitted(model.index)>0.5, 1, 0))))$CCR
ccr.app.kernel <- npconmode(bws=bw,newdata=train)$CCR.overall
## Conduct the M train/eval splits
for(m in 1:M) {
    ## Shuffle the data into independent training and evaluation samples.
    ii <- sample(1:n,replace=FALSE)
    train <- bwt[ii[1:n.train],]
    ## glm() can't deal with < all ftv cases
    while(length(unique(train$ftv))<length(unique(bwt$ftv))) {
        ii <- sample(1:n,replace=FALSE)
        train <- bwt[ii[1:n.train],]
    }
    eval <- bwt[ii[(n.train+1):n],]
    ## Extract the correct classification ratios for the independent
    ## evaluation data where we know the outcomes (update() refits the
    ## Logit model on train, the nonparametric model will
    ## automatically update taking train from the environment when it
    ## is called by predict()).
    ccr.linear[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linear),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.linint[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linint),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.BIC[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.BIC),
                                                type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.index[m] <- with(eval,CM(table(low,ifelse(predict(model.index,newdata=eval)>0.5,1,0))))$CCR
    ccr.kernel[m] <- npconmode(bws=bw,newdata=eval)$CCR.overall
}
## Conduct a paired t-test that the mean expected true CCR for each model 
## is equal versus the alternative that the kernel-based model has a significantly larger
## expected true CCR than the BIC-optimal Logit model.
p.linear <- t.test(ccr.kernel,ccr.linear,alternative="greater",paired=TRUE)$p.value
p.linint <- t.test(ccr.kernel,ccr.linint,alternative="greater",paired=TRUE)$p.value
p.BIC <- t.test(ccr.kernel,ccr.BIC,alternative="greater",paired=TRUE)$p.value
p.index <- t.test(ccr.kernel,ccr.index,alternative="greater",paired=TRUE)$p.value
p <- c(p.linear,p.linint,p.BIC,p.index,NA)
## Apparent performance
apparent <- c(ccr.app.linear,ccr.app.linint,ccr.app.BIC,ccr.app.index,ccr.app.kernel)
## Expected true performance
true <- c(mean(ccr.linear),mean(ccr.linint),mean(ccr.BIC),mean(ccr.index),mean(ccr.kernel))



foo <- data.frame(apparent,true,rank(-true),p)
rownames(foo) <- c("Par-Linear","Par-Linear-Int","Par-BIC","Semipar-KS","Nonpar")
colnames(foo) <- c("Apparent","Expected","Rank","$P$-value")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (higher values are preferred) and $P$-values from a test for equality of expected true performance. The cases considered are the kernel versus linear models, kernel versus linear with interaction models, kernel versus BIC-optimal models, and kernel versus the Klein and Spady models. Rejection of the null implies that the kernel-based model has significantly higher mean CCR on independent data.",
             digits=3,
             escape=FALSE)



## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(crs)
require(np)
data(wage1)
options(crs.messages=FALSE,np.messages=FALSE)
set.seed(42)
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), and number of train/eval splits (M).
n.eval <- 10
n <- nrow(wage1)
n.train <- n-n.eval
M <- 1000
## Create storage vectors.
aspe.linear <- numeric(M)
aspe.quadint <- numeric(M)
aspe.mw <- numeric(M)
aspe.vc <- numeric(M)
aspe.pl <- numeric(M)
aspe.si <- numeric(M)
aspe.lc <- numeric(M)
aspe.ll <- numeric(M)
aspe.glp <- numeric(M)
## Copy the full sample into the object train.
train <- wage1
## Fit the parametric, semiparametric, and nonparametric models for 
## the full sample.
##
## Linear 
lm.linear <- lm(lwage~female+
                married+
                educ+
                exper+
                tenure,
                data=train)
## Linear with interactions and quadratic in experience (common
## specification).
lm.quadint <- lm(lwage~(female+
                        married+
                        educ+
                        exper+
                        I(exper^2)+
                        tenure)^2,
                data=train)
## Murphy-Welch quartic specification                
lm.mw <- lm(lwage~female+
            married+
            educ+
            exper+
            I(exper^2)+
            I(exper^3)+
            I(exper^4)+
            tenure,
            data=train)
## Varying coefficient
model.vc <- npscoef(lwage~educ+
                    tenure+
                    exper+
                    expersq|female+married,
                    data=train)
## Partially linear
model.pl <- npplreg(lwage~female+
                    married+
                    educ+
                    tenure|exper,
                    data=train)
## Single-index (Ichimura)
model.si <- npindex(lwage~female+
                    married+
                    educ+
                    tenure+
                    exper+
                    expersq,
                    method="ichimura",
                    data=train)
## Local constant
bw.lc <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="lc",
                 bwmethod="cv.aic",
                 data=train)
model.lc <- npreg(bws=bw.lc)
## Local linear
bw.ll <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="ll",
                 bwmethod="cv.aic",
                 data=train)
model.ll <- npreg(bws=bw.ll)
## Generalized local polynomial
ghat <- npglpreg(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 cv.func="cv.aic",
                 data=train)
model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
## Apparent (in-sample) performance
aspe.app.linear <- with(train,mean((lwage-predict(lm.linear,newdata=train))^2))
aspe.app.quadint <- with(train,mean((lwage-predict(lm.quadint,newdata=train))^2))
aspe.app.mw <- with(train,mean((lwage-predict(lm.mw,newdata=train))^2))
aspe.app.vc <- with(train,mean((lwage-predict(model.vc,newdata=train))^2))
aspe.app.pl <- with(train,mean((lwage-predict(model.pl,newdata=train))^2))
aspe.app.si <- with(train,mean((lwage-predict(model.si,newdata=train))^2))
aspe.app.lc <- with(train,mean((lwage-predict(model.lc,newdata=train))^2))
aspe.app.ll <- with(train,mean((lwage-predict(model.ll,newdata=train))^2))
aspe.app.glp <- with(train,mean((lwage-predict(model.glp,newdata=train))^2))
## Conduct the M train/eval splits
for(m in 1:M) {
    ## We set the seed here to guarantee that the shuffles generated
    ## here and those in the Practitioner's Corner in Chapter 6 are
    ## identical.
    set.seed(m)
    ## Shuffle the data into independent training and evaluation
    ## samples.
    ii <- sample(1:n,replace=FALSE)
    train <- wage1[ii[1:n.train],]
    eval <- wage1[ii[(n.train+1):n],]
    ## Extract the APSEs for the independent evaluation data where we
    ## know the outcomes (update() refits the lm model on train).
    aspe.linear[m] <- with(eval,mean((lwage-predict(update(lm.linear),newdata=eval))^2))
    aspe.quadint[m] <- with(eval,mean((lwage-predict(update(lm.quadint),newdata=eval))^2))
    aspe.mw[m] <- with(eval,mean((lwage-predict(update(lm.mw),newdata=eval))^2))
    ## Calling the semi- and nonparametric functions with the existing
    ## bandwidth object re-estimates the model on the updated training data.
    model.vc.boot <- npscoef(bws=model.vc$bws)
    aspe.vc[m] <- with(eval,mean((lwage-predict(model.vc.boot,newdata=eval))^2))
    model.pl.boot <- npplreg(bws=model.pl$bw)
    aspe.pl[m] <- with(eval,mean((lwage-predict(model.pl.boot,newdata=eval))^2))
    model.si.boot <- npindex(bws=model.si$bws)
    aspe.si[m] <- with(eval,mean((lwage-predict(model.si.boot,newdata=eval))^2))
    model.lc <- npreg(bws=bw.lc)        
    aspe.lc[m] <- with(eval,mean((lwage-predict(model.lc,newdata=eval))^2))
    model.ll <- npreg(bws=bw.ll)
    aspe.ll[m] <- with(eval,mean((lwage-predict(model.ll,newdata=eval))^2))
    model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
    aspe.glp[m] <- with(eval,mean((lwage-predict(model.glp,newdata=eval))^2))
}
## Conduct a paired t-test that the mean expected true ASPE for each
## model is equal versus the alternative that the kernel has a
## significantly lower expected true ASPE than the MW-optimal lm
## model.
## LC versus parametric and semiparametric.
p.linear <- t.test(aspe.lc,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.lc,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.lc,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.lc,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.lc,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.lc,aspe.si,alternative="less",paired=TRUE)$p.value
p.lc <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## LL versus parametric and semiparametric.
p.linear <- t.test(aspe.ll,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.ll,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.ll,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.ll,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.ll,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.ll,aspe.si,alternative="less",paired=TRUE)$p.value
p.ll <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## GLP versus parametric and semiparametric.
p.linear <- t.test(aspe.glp,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.glp,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.glp,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.glp,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.glp,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.glp,aspe.si,alternative="less",paired=TRUE)$p.value
p.glp <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## Apparent performance.
apparent <- c(aspe.app.linear,aspe.app.quadint,aspe.app.mw,
              aspe.app.vc,aspe.app.pl,aspe.app.si,
              aspe.app.glp,aspe.app.ll,aspe.app.lc)
## Expected true performance.
true <- c(mean(aspe.linear),mean(aspe.quadint),mean(aspe.mw),
          mean(aspe.vc),mean(aspe.pl),mean(aspe.si),
          mean(aspe.glp),mean(aspe.ll),mean(aspe.lc))



foo <- data.frame(apparent,true,rank(true),p.glp,p.ll,p.lc)
rownames(foo) <- c("Par-Linear","Par-Quad-Int","Par-Murphy-Welch",
                   "Semipar-VC","Semipar-PL","Semipar-SI",
                   "Nonpar-GLP","Nonpar-LL","Nonpar-LC")
colnames(foo) <- c("Apparent","Expected","Rank",
                   "$P$-GLP","$P$-LL","$P$-LC")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (lower values are preferred) and $P$-values from a test for equality of expected true performance for the kernel versus the parametric and semiparametric models (rejection of the null implies the kernel model has significantly lower mean ASPE on independent data).",
             digits=4,
             escape=FALSE)


#| echo: true
#| eval: false
## ## Generate some data: sex (unordered categorical), income (ordered categorical),
## ## and height (numeric)
## n <- 100
## sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
## income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
## height <- rnorm(n,mean=150,sd=20)
## ## Note - by default these variables may not be the data types we desire
## class(sex);class(income);class(height)
## ## income is already numeric(), but sex and height are character()
## ## sex is categorical and unordered, so cast as type factor()
## sex <- factor(sex)
## ## income is categorical and ordered, but we need to ensure intended order (it
## ## will assume alphabetical ordering otherwise). Suppose you ignore it - let's
## ## see what happens when we just cast as type ordered() using defaults
## income <- ordered(income)
## ## The levels are in alphabetical order, which we don't want
## levels(income)
## ## We shall reorder the ordered factor levels as intended using levels=...
## income <- ordered(income,levels=c("Low","Middle","High"))
## levels(income)
## ## Check data types again
## class(sex);class(income);class(height)
## ## Note that with integers the default ordered() works fine
## x <- sample(c(2,5,4,3,1),n,replace=TRUE)
## x <- ordered(x)
## levels(x)


#| eval: true
#| echo: true
## Let's simulate a random numeric sample from the normal distribution
set.seed(42)
n <- 1000
x <- rnorm(n)
## Let's sort the data so we can graph x versus dnorm(x,...) using lines (type="l")
x <- sort(x)
## Since we simulated the data, let's plot the true, known, parametric density
## (we can't do this with actual data because the density of such data is, in
## general, unknown)
plot(x,dnorm(x,mean=mean(x),sd=sd(x)),type="l",ylab="Parametric Density Estimate",xlab="X")



pander::pander(shapiro.test(x))



data(faithful)
with(faithful,pander::pander(shapiro.test(eruptions)))



plot(density(faithful$eruptions),main="")
rug(faithful$eruptions)



plot(density(faithful$eruptions),main="")
with(faithful,lines(density(eruptions)$x,dnorm(density(eruptions)$x,mean=mean(eruptions),sd=sd(eruptions)),col=2,lty=2))
rug(faithful$eruptions)
legend("topleft",c("Nonparametric","Parametric"),lty=c(1,2),col=c(1,2),bty="n")



hist(faithful$eruptions,prob=TRUE,main="",xlab="Eruptions",breaks=20,xlim=c(1.25,5.5))
with(faithful,lines(density(eruptions)$x,fitted(npudens(tdat=eruptions,edat=density(eruptions)$x))))
rug(faithful$eruptions)


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions,data=faithful)
plot(fhat,neval=250,plot.errors.method="bootstrap")


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions+waiting,data=faithful)
plot(fhat,theta=330,xtrim=-0.05,neval=75,view="fixed",main="")


#| echo: true
n <- 250
set.seed(42)
sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
sex <- factor(sex)
phat <- npudens(~sex)
plot(phat,plot.errors.method="bootstrap")


#| echo: true
n <- 250
set.seed(42)
income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
income <- ordered(income,levels=c("Low","Middle","High"))
phat <- npudens(~income)
plot(phat,plot.errors.method="bootstrap")


#| label: tbl-wage1mixedtable
#| tbl-cap: Counts of number of dependants present in 526 households by cell
library(np)
library(plot3D)
data(wage1)
knitr::kable(with(wage1,t(data.frame(numdep=sort(unique(numdep)),counts=as.numeric(table(numdep))))),
             booktabs=TRUE,
             linesep="")


#| echo: true
#| eval: false
## library(np)
## library(plot3D)
## data(wage1)
## numdep.seq <- with(wage1,sort(unique(numdep)))
## lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
## wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
## bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
## fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## ## Hack since scatter3D converts ordered 0-6 to numeric 1-7
## scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
##           ylab="Log-Wage",
##           xlab="Number of Dependants",
##           zlab="Joint Density",
##           ticktype="detailed",
##           angle=15,
##           box=TRUE,
##           type="h",
##           grid=TRUE,
##           col="blue",
##           colkey=FALSE)


#| label: fig-wage1mixeddensity
#| fig-cap: Mixed-data bivariate kernel density estimate for the joint PDF of lwage (numeric) and numdeps (ordered)
data(wage1)
bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
numdep.seq <- with(wage1,sort(unique(numdep)))
lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## Hack since scatter3D converts ordered 0-6 to numeric 1-7
scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
          ylab="Log-Wage",
          xlab="Number of Dependants",
          zlab="Joint Density",
          ticktype="detailed",
          angle=15,
          box=TRUE,
          type="h",
          grid=TRUE,
          col="blue",
          colkey=FALSE)


rm(list=ls())
library(np)
library(car)
Prestige <- Prestige[order(Prestige$income),]
attach(Prestige)
## Plug-in
library(KernSmooth)
h <- dpill(income, prestige)
model.plugin <- locpoly(income, prestige, bandwidth = h)
## Under
bw <- npregbw(prestige~income,bws=c(sd(income)*length(income)^{-1/5}/10),regtype="ll",bandwidth.compute=FALSE)
fit.under <- fitted(npreg(bws=bw))
## Over
bw <- npregbw(prestige~income,bws=c(1000*sd(income)*length(income)^{-1/5}),regtype="ll",bandwidth.compute=FALSE)
fit.over <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.ls",regtype="ll")
fit.ls <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.aic",regtype="ll")
fit.aic <- fitted(npreg(bws=bw))
par(mfrow=c(2,2))
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Undersmoothed")
lines(income,fit.under)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Oversmoothed")
lines(income,fit.over)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Plug-In")
lines(model.plugin$x,model.plugin$y)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="AICc & CV-LS")
lines(income,fit.ls,lty=1)
lines(income,fit.aic,lty=2)


#| echo: true
library(np)
set.seed(42)
n <- 1000
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + rnorm(n,sd=0.25*sd(dgp))
ghat <- npreg(y~x)
plot(x,y,cex=0.5,col="grey",xlab="X",ylab="Y")
lines(x,dgp)
lines(x,fitted(ghat),col=2)
abline(ghat.ols <- lm(y~x),col=3)
legend("top",c("DGP","Kernel","OLS"),col=1:3,lty=1,bty="n")



pander::pander(summary(ghat.ols))


#| echo: true
plot(ghat,gradients=TRUE,neval=250)
lines(x,-2*pi*sin(2*pi*x),col=2)
abline(h=coef(ghat.ols)[2],col=3)
legend("topleft",c("Kernel ME","DGP ME","Linear ME"),col=1:3,lty=1,bty="n")


#| echo: true
library(np)
data(wage1)
## For details see ?wage1
ghat <- npreg(lwage ~ female + married + educ + exper + tenure, data=wage1, regtype="ll")
summary(ghat)


#| echo: true
## We run out of graph axis dimensions with > 2 predictors, so it is common to
## construct partial plots that plot the fitted model versus each predictor
## separately holding the off-axis predictors at, say, their median value (you
## can change this - see ?npplot and the argument xq)
par(mfrow=c(2,3))
plot(ghat,plot.errors.method="bootstrap")


#| echo: true
par(mfrow=c(2,3))
plot(ghat,gradients=TRUE,plot.errors.method="bootstrap")


#| echo: true
ghat.ols <- lm(lwage ~ female + married + educ + exper + tenure, data=wage1)
summary(ghat.ols)


#| echo: true
npsigtest(ghat)


#| echo: true
attach(wage1)
df <- data.frame(female = factor("Male", levels=levels(female)),
                 married = factor("Notmarried", levels=levels(married)),
                 educ = median(educ),
                 tenure = median(tenure),
                 exper = median(exper))
head(df)
predict(ghat, newdata=df)
## Or you could use ghat <- npreg(...,newdata=df) and fitted(ghat) 
## ghat <- npreg(lwage ~ female + married + educ + exper + tenure, 
##               data=wage1,
##               regtype="ll", 
##               newdata=df)
## fitted(ghat)
## Note - If `df` had multiple rows we would get a prediction for each
## row (i.e., a vector of predictions)

#| echo: true
#| eval: false
## ## Example 1 - create a simple time series that begins on the 2nd Quarter of
## ## 1959
## x <- ts(1:10, frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)
## ## Example 2 - create simulated data from an ARIMA model that begins on the 2nd
## ## Quarter of 1959
## x <- ts(arima.sim(n = 10, list(ar = c(0.8897, -0.4858),
##                                ma = c(-0.2279, 0.2488)),
##                   sd = sqrt(0.1796)),
##         frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)


#| echo: true
#| eval: false
## ## ldeaths is monthly deaths from bronchitis, emphysema and asthma in the UK,
## ## 1974–1979, both sexes
## ldeaths
## ## lag(x, k=1, ...), note the sign of k: a series lagged by a positive k, i.e.,
## ## y_{t+k}
## lag(ldeaths, 0) # returns original series y_t
## lag(ldeaths, -1) # returns series lagged once y_{t-1}
## lag(ldeaths, -2) # returns series lagged twice y_{t-2}
## ## In matrix form compare to cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))
## cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))


#| echo: true
#| eval: false
## diff(ldeaths)
## ## Manually compute first difference
## ldeaths[2:length(ldeaths)]-ldeaths[1:(length(ldeaths)-1)]


#| echo: true
#| eval: false
## library(np)
## ## Simulate a simple AR(1) process with lag coefficient 0.9
## set.seed(42)
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Using lag() in lm() calls will fail (same result regardless of lag, R^2=1,
## ## lag() is totally ignored so you are regressing y_t on y_t though naturally
## ## you may not suspect that is what would happen in this case)
## ghat.ols <- lm(y~lag(y,-1))
## summary(ghat.ols)
## ## But if you manually regress numeric y_t on y_{t-1} you get the expected
## ## result (try also ar.ols(y, order.max = 1) or arima(y, order = c(1,0,0))).
## ## Here we regress (y_t,...,y_2) on (y_{t-1},...,y_1) using numeric vectors
## ghat.ols <- lm(y[2:length(y)]~y[1:(length(y)-1)])
## summary(ghat.ols)
## ## npreg() does support ts() objects. In the example below we use the npreg()
## ## function to conduct simple OLS by noting that local linear regression with a
## ## large bandwidth is linear least squares, so the gradient will be a vector of
## ## constants equal to the least squares coefficient for lag(y,-1) - think of the
## ## lag() function and notation y_{t-1} or "why tee minus one" hence the -1 here
## ghat <- npreg(y~lag(y,-1),regtype="ll",bws=10^5,bandwidth.compute=FALSE,gradients=TRUE)
## tail(gradients(ghat))


#| echo: true
#| eval: false
## library(np)
## ## Simulate data from a stationary univariate AR(1) process
## set.seed(42)
## n <- 100
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Conduct local constant estimation of y on four lags
## ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## ## Re-fit the model using only 1 lag
## ghat.1.lag <- npreg(y~lag(y,-1))
## ## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## ## Plot the data and nonparametric fit
## plot(y,lwd=2)
## lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
## points(ts(fitted(ghat.1.lag),start=2),col=3)
## legend("topright",c("Series",
##                     "npreg() including 4 lags",
##                     "npreg() including 1 lag"),
##                     lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")


#| echo: false
#| label: fig-arimatsnpreg
#| fig-cap: Local constant estimation of a stationary univariate AR(1) time series, $Y_t=0.9Y_{t-1}+\epsilon_t$
library(np)
## Simulate data from a stationary univariate AR(1) process
set.seed(42)
n <- 100
y <- arima.sim(n = n, list(ar = c(0.9)))
## Conduct local constant estimation of y on four lags
ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## Re-fit the model using only 1 lag
ghat.1.lag <- npreg(y~lag(y,-1))
## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## Plot the data and nonparametric fit
plot(y,lwd=2)
lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
points(ts(fitted(ghat.1.lag),start=2),col=3)
legend("topright",c("Series",
                    "npreg() including 4 lags",
                    "npreg() including 1 lag"),
                    lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")



## Examine the model summary
summary(ghat.4.lags)


#| echo: true
## Conduct a nonparametric significance test
npsigtest(ghat.4.lags)


#| echo: true
## npglpreg() does not support time series objects, so we manually
## create the lagged series
library(crs)
options(crs.messages=FALSE)
y.lag.1 <- c(rep(NA,1),y[1:(n-1)])
y.lag.2 <- c(rep(NA,2),y[1:(n-2)])
y.lag.3 <- c(rep(NA,3),y[1:(n-3)])
y.lag.4 <- c(rep(NA,4),y[1:(n-4)])
model.glp <- npglpreg(y~y.lag.1+y.lag.2+y.lag.3+y.lag.4,nmulti=5,degree.max=5)


#| echo: true
summary(model.glp)


#| echo: false
model.lm <- lm(y~y.lag.1,subset=5:n)


#| echo: true
summary(model.lm)


#| echo: true
## Examine the first few fitted values from the nonparametric and
## linear AR(1) models
cbind(fitted(model.glp),fitted(model.lm))[1:10,]


#| echo: true
#| eval: false
## ## Fan & Yao's (1998) estimator, trim negative values of they occur
## ## Chen, Cheng & Peng's (2009) estimator is always positive but my
## ## Monte Carlo simulations show it is less efficient for this DGP
## library(np)
## set.seed(42)
## n <- 1000
## x <- sort(runif(n))
## sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
## dgp <- sin(2*pi*x)
## y <- dgp + rnorm(n,sd=sqrt(sigma.var))
## model <- npreg(y~x,regtype="ll")
## r <- residuals(model)**2
## ## Fan and Yao's (1998) estimator with trimming if needed
## var.fy <- fitted(npreg(r~x,regtype="ll"))
## var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## ## Chen, Cheng and Peng's (2009) estimator
## log.r <- log(r+1/n) ## Avoids log(0)
## V.hat <- fitted(npreg(log.r~x,regtype="ll"))
## d.hat <- 1/mean(r*exp(-V.hat))
## var.ccp <- exp(V.hat)/d.hat
## par(mfrow=c(1,2))
## plot(x,y,cex=.25,col="grey")
## lines(x,dgp,col=1,lty=1,lwd=2)
## lines(x,fitted(model),col=2,lty=2,lwd=2)
## legend("topleft",c("g(x)=sin(2 pi x)",
##                    "Nonparametric Estimate"),
##        col=1:2,
##        lty=1:2,
##        lwd=c(2,2),
##        bty="n")
## ylim=c(min(r),quantile(r,0.95))
## plot(x,r,ylim=ylim,cex=.25,col="grey")
## lines(x,sigma.var,col=1,lty=1,lwd=2)
## lines(x,var.fy,col=2,lty=2,lwd=2)
## lines(x,var.ccp,col=3,lty=4,lwd=2)
## legend("topleft",c("volatility=(2 pi x)^2/10",
##                    "Fan and Yao (1998)",
##                    "Chen et al. (2009)"),
##        col=1:3,
##        lty=c(1,2,4),
##        lwd=c(2,2,2),
##        bty="n")


#| label: fig-fanyaosim
#| fig-cap: Fan and Yao's (1998) and Chen et al.'s (2009) conditional variance estimators
#| echo: false
par(mfrow=c(1,2))
## Fan & Yao's (1998) estimator, trim negative values of they occur
## Chen, Cheng & Peng's (2009) estimator is always positive but my
## Monte Carlo simulations show it is less efficient for this DGP
set.seed(42)
n <- 1000
x <- sort(runif(n))
sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
dgp <- sin(2*pi*x)
y <- dgp + rnorm(n,sd=sqrt(sigma.var))
model <- npreg(y~x,regtype="ll")
r <- residuals(model)**2
## Fan and Yao's (1998) estimator with trimming if needed
var.fy <- fitted(npreg(r~x,regtype="ll"))
var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## Chen, Cheng and Peng's (2009) estimator
log.r <- log(r+1/n) ## Avoids log(0)
V.hat <- fitted(npreg(log.r~x,regtype="ll"))
d.hat <- 1/mean(r*exp(-V.hat))
var.ccp <- exp(V.hat)/d.hat
plot(x,y,cex=.25,col="grey")
lines(x,dgp,col=1,lty=1,lwd=2)
lines(x,fitted(model),col=2,lty=2,lwd=2)
legend("topleft",c("g(x)=sin(2 pi x)", 
                   "Nonparametric Estimate"),
       col=1:2,
       lty=1:2,
       lwd=c(2,2),
       bty="n")
ylim=c(min(r),quantile(r,0.95))
plot(x,r,ylim=ylim,cex=.25,col="grey")
lines(x,sigma.var,col=1,lty=1,lwd=2)
lines(x,var.fy,col=2,lty=2,lwd=2)
lines(x,var.ccp,col=3,lty=4,lwd=2)
legend("topleft",c("volatility=(2 pi x)^2/10",
                   "Fan and Yao (1998)",
                   "Chen et al. (2009)"),
       col=1:3,
       lty=c(1,2,4),
       lwd=c(2,2,2),
       bty="n")


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## ## Compute the pooled first-step estimator
## model.pooled <- npreg(y~x,ckerorder=4)
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~x,ckerorder=2)


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## model.fe <- npreg(y~x+factor(id))


#| echo: true
#| eval: false
## library(Ecdat)
## data(Airline)
## ?Airline
## library(np)
## options(np.messages=FALSE)
## attach(Airline)
## 
## ## Dependent variable is (log) cost, predictors (log) output, fuel price, load factor, year
## lcost <- as.numeric(log(cost))
## loutput <- as.numeric(log(output))
## lpf <- as.numeric(log(pf))
## lf <- as.numeric(lf)
## year <- ordered(year)
## airline <- factor(airline)
## 
## ## Airline specific fixed effects
## 
## model.fe <- npreg(lcost~loutput+lpf+lf+year+airline,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine",
##                   gradients=TRUE)
## 
## summary(model.fe)
## summary(model.fe)
## ## Plot partial means
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)
## 
## ## Random effects
## 
## ## Compute the pooled first-step estimator
## 
## id <- airline
## n <- length(unique(airline))
## T <- length(unique(year))
## model.pooled <- npreg(lcost~loutput+lpf+lf+year,
##                       ckerorder=4,
##                       regtype="ll",
##                       bwmethod="cv.aic",
##                       ukertype="liracine",
##                       okertype="liracine")
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~loutput+lpf+lf+year,
##                   ckerorder=2,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine")
## 
## summary(model.re)
## ## Plot partial means
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)

#| echo: true
#| eval: true
library(np)
data(wage1)
model.lm <- lm(lwage~female+married+educ+tenure+exper+I(exper^2),data=wage1)
knitr::kable(coef(model.lm),col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.pl <- npplreg(lwage~female+married+educ+tenure|exper,
                    data=wage1)
knitr::kable(coef(model.pl),digits=7,col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.scoef <- npscoef(lwage~educ+tenure+exper+expersq|female+married,
                       data=wage1,
                       betas=TRUE)
summary(model.scoef)


#| echo: true
#| eval: true
## The smooth coefficients are vectors, one for each predictor in X (they "vary"
## with Z), so we compute the means of the columns of these coefficients, one for
## each X predictor
knitr::kable(colMeans(coef(model.scoef)),digits=7,col.names="Avg. Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.index <- npindex(lwage~educ+
                       tenure+
                       exper+
                       expersq+
                       female+
                       married,
                       method="ichimura",
                       data=wage1)
summary(model.index)


#| echo: true
#| eval: true
library(MASS)
data("birthwt")
model.logit <- glm(low~factor(smoke)+
                   factor(race)+
                   factor(ht)+
                   factor(ui)+
                   ordered(ftv)+
                   age+
                   lwt,
                   family=binomial(link=logit),data=birthwt)
cm.logit <- with(birthwt,table(low, ifelse(fitted(model.logit)>0.5, 1, 0)))
ccr.logit <- sum(diag(cm.logit))/sum(cm.logit)
knitr::kable(cm.logit)


#| echo: true
#| eval: true
library(MASS)
data(birthwt)
library(np)
model.index <- npindex(low~factor(smoke)+
                       factor(race)+
                       factor(ht)+
                       factor(ui)+
                       ordered(ftv)+
                       age+
                       lwt,
                       method="kleinspady",
                       data=birthwt)
cm.index <- with(birthwt,table(low, ifelse(fitted(model.index)>0.5, 1, 0)))
ccr.index <- sum(diag(cm.index))/sum(cm.index)
knitr::kable(cm.index)


require(np)
require(MASS)
## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(np)
options(np.messages=FALSE)
set.seed(42)
## This function takes a confusion matrix and formats it correctly if
## it is unbalanced and returns the CCR as well.
CM <- function(cm) {
  factor.values.eval <- colnames(cm)
  CM <- matrix(0,nrow(cm),nrow(cm))
  rownames(CM) <- rownames(cm)
  colnames(CM) <- rownames(cm)
  for(i in 1:ncol(cm)) CM[,(1:nrow(cm))[rownames(cm)==factor.values.eval[i]]] <- cm[,i]
  return(list(CM=CM,CCR=sum(diag(CM))/sum(CM)))
}
## Load the birthwt data
library(MASS)
data(birthwt)
## Create a data frame that has up to 4th order polynomials. We will
## use the BIC-optimal parametric model that will allow for
## interactions and some nonlinearity.
bwt <- with(birthwt,data.frame(low=factor(low),
                               race=factor(race),
                               smoke=factor(smoke),
                               ht=factor(ht),
                               ui=factor(ui),
                               ftv=ordered(ftv),
                               age=age,
                               agesq=age**2,
                               agecu=age**3,
                               agequ=age**4,
                               lwt=lwt,
                               lwtsq=lwt**2,
                               lwtcu=lwt**3,
                               lwtqu=lwt**4))
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), number of multistarts for bandwidth selection (nmulti)
## and number of train/eval splits (M)
n.eval <- 10
n <- nrow(bwt)
n.train <- n-n.eval
M <- 1000
## Create storage vectors
ccr.linear <- numeric(M)
ccr.linint <- numeric(M)
ccr.BIC <- numeric(M)
ccr.kernel <- numeric(M)
ccr.index <- numeric(M)
## Copy the full sample into the object train
train <- bwt
## Fit the parametric Logit model for the full sample.
##
## Linear 
logit.linear <- glm(low~
                    smoke+
                    race+
                    ht+
                    ui+
                    age+
                    lwt+
                    ftv,
                    family=binomial(link=logit),
                    data=train)
## Linear Logit model with interactions
logit.linint <- glm(low~
                    (smoke+
                     race+
                     ht+
                     ui+
                     age+
                     lwt+
                     ftv)^2,
                    family=binomial(link=logit),
                    data=train)
## BIC-optimal Logit model
logit.BIC <- glm(low ~ .,
                 family = binomial(link=logit),
                 data = train)
logit.BIC <- stepAIC(logit.BIC, ~ .^3,
                     trace=FALSE,
                     k=log(nrow(birthwt)))
## Klein spady estimator
model.index <- npindex(low~smoke+
                       race+
                       ht+
                       ui+
                       ftv+
                       age+
                       lwt,
                       method="kleinspady",
                       data=train)
## Get the bandwidths for the nonparametric model for the full sample.
bw <- npcdensbw(low~
                smoke+
                race+
                ht+
                ui+
                age+
                lwt+
                ftv,
                data=train)
## Apparent (in-sample) performance
ccr.app.linear <- with(train,CM(table(low,ifelse(predict(update(logit.linear),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.linint <- with(train,CM(table(low,ifelse(predict(update(logit.linint),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.BIC <- with(train,CM(table(low,ifelse(predict(update(logit.BIC),
                                              type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.index <- with(train,CM(table(low,ifelse(fitted(model.index)>0.5, 1, 0))))$CCR
ccr.app.kernel <- npconmode(bws=bw,newdata=train)$CCR.overall
## Conduct the M train/eval splits
for(m in 1:M) {
    ## Shuffle the data into independent training and evaluation samples.
    ii <- sample(1:n,replace=FALSE)
    train <- bwt[ii[1:n.train],]
    ## glm() can't deal with < all ftv cases
    while(length(unique(train$ftv))<length(unique(bwt$ftv))) {
        ii <- sample(1:n,replace=FALSE)
        train <- bwt[ii[1:n.train],]
    }
    eval <- bwt[ii[(n.train+1):n],]
    ## Extract the correct classification ratios for the independent
    ## evaluation data where we know the outcomes (update() refits the
    ## Logit model on train, the nonparametric model will
    ## automatically update taking train from the environment when it
    ## is called by predict()).
    ccr.linear[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linear),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.linint[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linint),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.BIC[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.BIC),
                                                type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.index[m] <- with(eval,CM(table(low,ifelse(predict(model.index,newdata=eval)>0.5,1,0))))$CCR
    ccr.kernel[m] <- npconmode(bws=bw,newdata=eval)$CCR.overall
}
## Conduct a paired t-test that the mean expected true CCR for each model 
## is equal versus the alternative that the kernel-based model has a significantly larger
## expected true CCR than the BIC-optimal Logit model.
p.linear <- t.test(ccr.kernel,ccr.linear,alternative="greater",paired=TRUE)$p.value
p.linint <- t.test(ccr.kernel,ccr.linint,alternative="greater",paired=TRUE)$p.value
p.BIC <- t.test(ccr.kernel,ccr.BIC,alternative="greater",paired=TRUE)$p.value
p.index <- t.test(ccr.kernel,ccr.index,alternative="greater",paired=TRUE)$p.value
p <- c(p.linear,p.linint,p.BIC,p.index,NA)
## Apparent performance
apparent <- c(ccr.app.linear,ccr.app.linint,ccr.app.BIC,ccr.app.index,ccr.app.kernel)
## Expected true performance
true <- c(mean(ccr.linear),mean(ccr.linint),mean(ccr.BIC),mean(ccr.index),mean(ccr.kernel))



foo <- data.frame(apparent,true,rank(-true),p)
rownames(foo) <- c("Par-Linear","Par-Linear-Int","Par-BIC","Semipar-KS","Nonpar")
colnames(foo) <- c("Apparent","Expected","Rank","$P$-value")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (higher values are preferred) and $P$-values from a test for equality of expected true performance. The cases considered are the kernel versus linear models, kernel versus linear with interaction models, kernel versus BIC-optimal models, and kernel versus the Klein and Spady models. Rejection of the null implies that the kernel-based model has significantly higher mean CCR on independent data.",
             digits=3,
             escape=FALSE)



## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(crs)
require(np)
data(wage1)
options(crs.messages=FALSE,np.messages=FALSE)
set.seed(42)
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), and number of train/eval splits (M).
n.eval <- 10
n <- nrow(wage1)
n.train <- n-n.eval
M <- 1000
## Create storage vectors.
aspe.linear <- numeric(M)
aspe.quadint <- numeric(M)
aspe.mw <- numeric(M)
aspe.vc <- numeric(M)
aspe.pl <- numeric(M)
aspe.si <- numeric(M)
aspe.lc <- numeric(M)
aspe.ll <- numeric(M)
aspe.glp <- numeric(M)
## Copy the full sample into the object train.
train <- wage1
## Fit the parametric, semiparametric, and nonparametric models for 
## the full sample.
##
## Linear 
lm.linear <- lm(lwage~female+
                married+
                educ+
                exper+
                tenure,
                data=train)
## Linear with interactions and quadratic in experience (common
## specification).
lm.quadint <- lm(lwage~(female+
                        married+
                        educ+
                        exper+
                        I(exper^2)+
                        tenure)^2,
                data=train)
## Murphy-Welch quartic specification                
lm.mw <- lm(lwage~female+
            married+
            educ+
            exper+
            I(exper^2)+
            I(exper^3)+
            I(exper^4)+
            tenure,
            data=train)
## Varying coefficient
model.vc <- npscoef(lwage~educ+
                    tenure+
                    exper+
                    expersq|female+married,
                    data=train)
## Partially linear
model.pl <- npplreg(lwage~female+
                    married+
                    educ+
                    tenure|exper,
                    data=train)
## Single-index (Ichimura)
model.si <- npindex(lwage~female+
                    married+
                    educ+
                    tenure+
                    exper+
                    expersq,
                    method="ichimura",
                    data=train)
## Local constant
bw.lc <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="lc",
                 bwmethod="cv.aic",
                 data=train)
model.lc <- npreg(bws=bw.lc)
## Local linear
bw.ll <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="ll",
                 bwmethod="cv.aic",
                 data=train)
model.ll <- npreg(bws=bw.ll)
## Generalized local polynomial
ghat <- npglpreg(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 cv.func="cv.aic",
                 data=train)
model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
## Apparent (in-sample) performance
aspe.app.linear <- with(train,mean((lwage-predict(lm.linear,newdata=train))^2))
aspe.app.quadint <- with(train,mean((lwage-predict(lm.quadint,newdata=train))^2))
aspe.app.mw <- with(train,mean((lwage-predict(lm.mw,newdata=train))^2))
aspe.app.vc <- with(train,mean((lwage-predict(model.vc,newdata=train))^2))
aspe.app.pl <- with(train,mean((lwage-predict(model.pl,newdata=train))^2))
aspe.app.si <- with(train,mean((lwage-predict(model.si,newdata=train))^2))
aspe.app.lc <- with(train,mean((lwage-predict(model.lc,newdata=train))^2))
aspe.app.ll <- with(train,mean((lwage-predict(model.ll,newdata=train))^2))
aspe.app.glp <- with(train,mean((lwage-predict(model.glp,newdata=train))^2))
## Conduct the M train/eval splits
for(m in 1:M) {
    ## We set the seed here to guarantee that the shuffles generated
    ## here and those in the Practitioner's Corner in Chapter 6 are
    ## identical.
    set.seed(m)
    ## Shuffle the data into independent training and evaluation
    ## samples.
    ii <- sample(1:n,replace=FALSE)
    train <- wage1[ii[1:n.train],]
    eval <- wage1[ii[(n.train+1):n],]
    ## Extract the APSEs for the independent evaluation data where we
    ## know the outcomes (update() refits the lm model on train).
    aspe.linear[m] <- with(eval,mean((lwage-predict(update(lm.linear),newdata=eval))^2))
    aspe.quadint[m] <- with(eval,mean((lwage-predict(update(lm.quadint),newdata=eval))^2))
    aspe.mw[m] <- with(eval,mean((lwage-predict(update(lm.mw),newdata=eval))^2))
    ## Calling the semi- and nonparametric functions with the existing
    ## bandwidth object re-estimates the model on the updated training data.
    model.vc.boot <- npscoef(bws=model.vc$bws)
    aspe.vc[m] <- with(eval,mean((lwage-predict(model.vc.boot,newdata=eval))^2))
    model.pl.boot <- npplreg(bws=model.pl$bw)
    aspe.pl[m] <- with(eval,mean((lwage-predict(model.pl.boot,newdata=eval))^2))
    model.si.boot <- npindex(bws=model.si$bws)
    aspe.si[m] <- with(eval,mean((lwage-predict(model.si.boot,newdata=eval))^2))
    model.lc <- npreg(bws=bw.lc)        
    aspe.lc[m] <- with(eval,mean((lwage-predict(model.lc,newdata=eval))^2))
    model.ll <- npreg(bws=bw.ll)
    aspe.ll[m] <- with(eval,mean((lwage-predict(model.ll,newdata=eval))^2))
    model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
    aspe.glp[m] <- with(eval,mean((lwage-predict(model.glp,newdata=eval))^2))
}
## Conduct a paired t-test that the mean expected true ASPE for each
## model is equal versus the alternative that the kernel has a
## significantly lower expected true ASPE than the MW-optimal lm
## model.
## LC versus parametric and semiparametric.
p.linear <- t.test(aspe.lc,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.lc,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.lc,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.lc,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.lc,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.lc,aspe.si,alternative="less",paired=TRUE)$p.value
p.lc <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## LL versus parametric and semiparametric.
p.linear <- t.test(aspe.ll,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.ll,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.ll,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.ll,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.ll,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.ll,aspe.si,alternative="less",paired=TRUE)$p.value
p.ll <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## GLP versus parametric and semiparametric.
p.linear <- t.test(aspe.glp,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.glp,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.glp,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.glp,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.glp,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.glp,aspe.si,alternative="less",paired=TRUE)$p.value
p.glp <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## Apparent performance.
apparent <- c(aspe.app.linear,aspe.app.quadint,aspe.app.mw,
              aspe.app.vc,aspe.app.pl,aspe.app.si,
              aspe.app.glp,aspe.app.ll,aspe.app.lc)
## Expected true performance.
true <- c(mean(aspe.linear),mean(aspe.quadint),mean(aspe.mw),
          mean(aspe.vc),mean(aspe.pl),mean(aspe.si),
          mean(aspe.glp),mean(aspe.ll),mean(aspe.lc))



foo <- data.frame(apparent,true,rank(true),p.glp,p.ll,p.lc)
rownames(foo) <- c("Par-Linear","Par-Quad-Int","Par-Murphy-Welch",
                   "Semipar-VC","Semipar-PL","Semipar-SI",
                   "Nonpar-GLP","Nonpar-LL","Nonpar-LC")
colnames(foo) <- c("Apparent","Expected","Rank",
                   "$P$-GLP","$P$-LL","$P$-LC")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (lower values are preferred) and $P$-values from a test for equality of expected true performance for the kernel versus the parametric and semiparametric models (rejection of the null implies the kernel model has significantly lower mean ASPE on independent data).",
             digits=4,
             escape=FALSE)


#| echo: true
#| eval: false
## ## Generate some data: sex (unordered categorical), income (ordered categorical),
## ## and height (numeric)
## n <- 100
## sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
## income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
## height <- rnorm(n,mean=150,sd=20)
## ## Note - by default these variables may not be the data types we desire
## class(sex);class(income);class(height)
## ## income is already numeric(), but sex and height are character()
## ## sex is categorical and unordered, so cast as type factor()
## sex <- factor(sex)
## ## income is categorical and ordered, but we need to ensure intended order (it
## ## will assume alphabetical ordering otherwise). Suppose you ignore it - let's
## ## see what happens when we just cast as type ordered() using defaults
## income <- ordered(income)
## ## The levels are in alphabetical order, which we don't want
## levels(income)
## ## We shall reorder the ordered factor levels as intended using levels=...
## income <- ordered(income,levels=c("Low","Middle","High"))
## levels(income)
## ## Check data types again
## class(sex);class(income);class(height)
## ## Note that with integers the default ordered() works fine
## x <- sample(c(2,5,4,3,1),n,replace=TRUE)
## x <- ordered(x)
## levels(x)


#| eval: true
#| echo: true
## Let's simulate a random numeric sample from the normal distribution
set.seed(42)
n <- 1000
x <- rnorm(n)
## Let's sort the data so we can graph x versus dnorm(x,...) using lines (type="l")
x <- sort(x)
## Since we simulated the data, let's plot the true, known, parametric density
## (we can't do this with actual data because the density of such data is, in
## general, unknown)
plot(x,dnorm(x,mean=mean(x),sd=sd(x)),type="l",ylab="Parametric Density Estimate",xlab="X")



pander::pander(shapiro.test(x))



data(faithful)
with(faithful,pander::pander(shapiro.test(eruptions)))



plot(density(faithful$eruptions),main="")
rug(faithful$eruptions)



plot(density(faithful$eruptions),main="")
with(faithful,lines(density(eruptions)$x,dnorm(density(eruptions)$x,mean=mean(eruptions),sd=sd(eruptions)),col=2,lty=2))
rug(faithful$eruptions)
legend("topleft",c("Nonparametric","Parametric"),lty=c(1,2),col=c(1,2),bty="n")



hist(faithful$eruptions,prob=TRUE,main="",xlab="Eruptions",breaks=20,xlim=c(1.25,5.5))
with(faithful,lines(density(eruptions)$x,fitted(npudens(tdat=eruptions,edat=density(eruptions)$x))))
rug(faithful$eruptions)


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions,data=faithful)
plot(fhat,neval=250,plot.errors.method="bootstrap")


#| echo: true
library(np)
data(faithful)
fhat <- npudens(~eruptions+waiting,data=faithful)
plot(fhat,theta=330,xtrim=-0.05,neval=75,view="fixed",main="")


#| echo: true
n <- 250
set.seed(42)
sex <- sample(c("Female","Male"),n,replace=TRUE,prob=c(.4,.6))
sex <- factor(sex)
phat <- npudens(~sex)
plot(phat,plot.errors.method="bootstrap")


#| echo: true
n <- 250
set.seed(42)
income <- sample(c("Low","Middle","High"),n,replace=TRUE,prob=c(.3,.5,.2))
income <- ordered(income,levels=c("Low","Middle","High"))
phat <- npudens(~income)
plot(phat,plot.errors.method="bootstrap")


#| label: tbl-wage1mixedtable
#| tbl-cap: Counts of number of dependants present in 526 households by cell
library(np)
library(plot3D)
data(wage1)
knitr::kable(with(wage1,t(data.frame(numdep=sort(unique(numdep)),counts=as.numeric(table(numdep))))),
             booktabs=TRUE,
             linesep="")


#| echo: true
#| eval: false
## library(np)
## library(plot3D)
## data(wage1)
## numdep.seq <- with(wage1,sort(unique(numdep)))
## lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
## wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
## bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
## fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## ## Hack since scatter3D converts ordered 0-6 to numeric 1-7
## scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
##           ylab="Log-Wage",
##           xlab="Number of Dependants",
##           zlab="Joint Density",
##           ticktype="detailed",
##           angle=15,
##           box=TRUE,
##           type="h",
##           grid=TRUE,
##           col="blue",
##           colkey=FALSE)


#| label: fig-wage1mixeddensity
#| fig-cap: Mixed-data bivariate kernel density estimate for the joint PDF of lwage (numeric) and numdeps (ordered)
data(wage1)
bw <- npudensbw(~lwage+ordered(numdep),data=wage1)
numdep.seq <- with(wage1,sort(unique(numdep)))
lwage.seq <- with(wage1,seq(min(lwage),max(lwage),length=50))
wage1.eval <- expand.grid(numdep=ordered(numdep.seq),lwage=lwage.seq)
fhat <- fitted(npudens(bws=bw,newdata=wage1.eval))
## Hack since scatter3D converts ordered 0-6 to numeric 1-7
scatter3D(as.numeric(wage1.eval[,1])-1,wage1.eval[,2],fhat,
          ylab="Log-Wage",
          xlab="Number of Dependants",
          zlab="Joint Density",
          ticktype="detailed",
          angle=15,
          box=TRUE,
          type="h",
          grid=TRUE,
          col="blue",
          colkey=FALSE)


rm(list=ls())
library(np)
library(car)
Prestige <- Prestige[order(Prestige$income),]
attach(Prestige)
## Plug-in
library(KernSmooth)
h <- dpill(income, prestige)
model.plugin <- locpoly(income, prestige, bandwidth = h)
## Under
bw <- npregbw(prestige~income,bws=c(sd(income)*length(income)^{-1/5}/10),regtype="ll",bandwidth.compute=FALSE)
fit.under <- fitted(npreg(bws=bw))
## Over
bw <- npregbw(prestige~income,bws=c(1000*sd(income)*length(income)^{-1/5}),regtype="ll",bandwidth.compute=FALSE)
fit.over <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.ls",regtype="ll")
fit.ls <- fitted(npreg(bws=bw))
## LL CV
bw <- npregbw(prestige~income,bwmethod="cv.aic",regtype="ll")
fit.aic <- fitted(npreg(bws=bw))
par(mfrow=c(2,2))
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Undersmoothed")
lines(income,fit.under)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Oversmoothed")
lines(income,fit.over)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="Plug-In")
lines(model.plugin$x,model.plugin$y)
plot(income,prestige,xlab="Income (K)", ylab="Prestige",main="AICc & CV-LS")
lines(income,fit.ls,lty=1)
lines(income,fit.aic,lty=2)


#| echo: true
library(np)
set.seed(42)
n <- 1000
x <- sort(runif(n))
dgp <- cos(2*pi*x)
y <- dgp + rnorm(n,sd=0.25*sd(dgp))
ghat <- npreg(y~x)
plot(x,y,cex=0.5,col="grey",xlab="X",ylab="Y")
lines(x,dgp)
lines(x,fitted(ghat),col=2)
abline(ghat.ols <- lm(y~x),col=3)
legend("top",c("DGP","Kernel","OLS"),col=1:3,lty=1,bty="n")



pander::pander(summary(ghat.ols))


#| echo: true
plot(ghat,gradients=TRUE,neval=250)
lines(x,-2*pi*sin(2*pi*x),col=2)
abline(h=coef(ghat.ols)[2],col=3)
legend("topleft",c("Kernel ME","DGP ME","Linear ME"),col=1:3,lty=1,bty="n")


#| echo: true
library(np)
data(wage1)
## For details see ?wage1
ghat <- npreg(lwage ~ female + married + educ + exper + tenure, data=wage1, regtype="ll")
summary(ghat)


#| echo: true
## We run out of graph axis dimensions with > 2 predictors, so it is common to
## construct partial plots that plot the fitted model versus each predictor
## separately holding the off-axis predictors at, say, their median value (you
## can change this - see ?npplot and the argument xq)
par(mfrow=c(2,3))
plot(ghat,plot.errors.method="bootstrap")


#| echo: true
par(mfrow=c(2,3))
plot(ghat,gradients=TRUE,plot.errors.method="bootstrap")


#| echo: true
ghat.ols <- lm(lwage ~ female + married + educ + exper + tenure, data=wage1)
summary(ghat.ols)


#| echo: true
npsigtest(ghat)


#| echo: true
attach(wage1)
df <- data.frame(female = factor("Male", levels=levels(female)),
                 married = factor("Notmarried", levels=levels(married)),
                 educ = median(educ),
                 tenure = median(tenure),
                 exper = median(exper))
head(df)
predict(ghat, newdata=df)
## Or you could use ghat <- npreg(...,newdata=df) and fitted(ghat) 
## ghat <- npreg(lwage ~ female + married + educ + exper + tenure, 
##               data=wage1,
##               regtype="ll", 
##               newdata=df)
## fitted(ghat)
## Note - If `df` had multiple rows we would get a prediction for each
## row (i.e., a vector of predictions)

#| echo: true
#| eval: false
## ## Example 1 - create a simple time series that begins on the 2nd Quarter of
## ## 1959
## x <- ts(1:10, frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)
## ## Example 2 - create simulated data from an ARIMA model that begins on the 2nd
## ## Quarter of 1959
## x <- ts(arima.sim(n = 10, list(ar = c(0.8897, -0.4858),
##                                ma = c(-0.2279, 0.2488)),
##                   sd = sqrt(0.1796)),
##         frequency = 4, start = c(1959, 2))
## print(x)
## plot(x)


#| echo: true
#| eval: false
## ## ldeaths is monthly deaths from bronchitis, emphysema and asthma in the UK,
## ## 1974–1979, both sexes
## ldeaths
## ## lag(x, k=1, ...), note the sign of k: a series lagged by a positive k, i.e.,
## ## y_{t+k}
## lag(ldeaths, 0) # returns original series y_t
## lag(ldeaths, -1) # returns series lagged once y_{t-1}
## lag(ldeaths, -2) # returns series lagged twice y_{t-2}
## ## In matrix form compare to cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))
## cbind(ldeaths,lag(ldeaths,-1),lag(ldeaths,-2))


#| echo: true
#| eval: false
## diff(ldeaths)
## ## Manually compute first difference
## ldeaths[2:length(ldeaths)]-ldeaths[1:(length(ldeaths)-1)]


#| echo: true
#| eval: false
## library(np)
## ## Simulate a simple AR(1) process with lag coefficient 0.9
## set.seed(42)
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Using lag() in lm() calls will fail (same result regardless of lag, R^2=1,
## ## lag() is totally ignored so you are regressing y_t on y_t though naturally
## ## you may not suspect that is what would happen in this case)
## ghat.ols <- lm(y~lag(y,-1))
## summary(ghat.ols)
## ## But if you manually regress numeric y_t on y_{t-1} you get the expected
## ## result (try also ar.ols(y, order.max = 1) or arima(y, order = c(1,0,0))).
## ## Here we regress (y_t,...,y_2) on (y_{t-1},...,y_1) using numeric vectors
## ghat.ols <- lm(y[2:length(y)]~y[1:(length(y)-1)])
## summary(ghat.ols)
## ## npreg() does support ts() objects. In the example below we use the npreg()
## ## function to conduct simple OLS by noting that local linear regression with a
## ## large bandwidth is linear least squares, so the gradient will be a vector of
## ## constants equal to the least squares coefficient for lag(y,-1) - think of the
## ## lag() function and notation y_{t-1} or "why tee minus one" hence the -1 here
## ghat <- npreg(y~lag(y,-1),regtype="ll",bws=10^5,bandwidth.compute=FALSE,gradients=TRUE)
## tail(gradients(ghat))


#| echo: true
#| eval: false
## library(np)
## ## Simulate data from a stationary univariate AR(1) process
## set.seed(42)
## n <- 100
## y <- arima.sim(n = n, list(ar = c(0.9)))
## ## Conduct local constant estimation of y on four lags
## ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## ## Re-fit the model using only 1 lag
## ghat.1.lag <- npreg(y~lag(y,-1))
## ## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## ## Plot the data and nonparametric fit
## plot(y,lwd=2)
## lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
## points(ts(fitted(ghat.1.lag),start=2),col=3)
## legend("topright",c("Series",
##                     "npreg() including 4 lags",
##                     "npreg() including 1 lag"),
##                     lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")


#| echo: false
#| label: fig-arimatsnpreg
#| fig-cap: Local constant estimation of a stationary univariate AR(1) time series, $Y_t=0.9Y_{t-1}+\epsilon_t$
library(np)
## Simulate data from a stationary univariate AR(1) process
set.seed(42)
n <- 100
y <- arima.sim(n = n, list(ar = c(0.9)))
## Conduct local constant estimation of y on four lags
ghat.4.lags <- npreg(y~lag(y,-1)+lag(y,-2)+lag(y,-3)+lag(y,-4))
## Re-fit the model using only 1 lag
ghat.1.lag <- npreg(y~lag(y,-1))
## cbind(ts(fitted(ghat.4.lags),start=5),ts(fitted(ghat.1.lag),start=2))
## Plot the data and nonparametric fit
plot(y,lwd=2)
lines(ts(fitted(ghat.4.lags),start=5),col=2,lty=2,lwd=2)
points(ts(fitted(ghat.1.lag),start=2),col=3)
legend("topright",c("Series",
                    "npreg() including 4 lags",
                    "npreg() including 1 lag"),
                    lty=c(1,2,NA),col=1:3,lwd=c(2,2,NA),pch=c(NA,NA,1),bty="n")



## Examine the model summary
summary(ghat.4.lags)


#| echo: true
## Conduct a nonparametric significance test
npsigtest(ghat.4.lags)


#| echo: true
## npglpreg() does not support time series objects, so we manually
## create the lagged series
library(crs)
options(crs.messages=FALSE)
y.lag.1 <- c(rep(NA,1),y[1:(n-1)])
y.lag.2 <- c(rep(NA,2),y[1:(n-2)])
y.lag.3 <- c(rep(NA,3),y[1:(n-3)])
y.lag.4 <- c(rep(NA,4),y[1:(n-4)])
model.glp <- npglpreg(y~y.lag.1+y.lag.2+y.lag.3+y.lag.4,nmulti=5,degree.max=5)


#| echo: true
summary(model.glp)


#| echo: false
model.lm <- lm(y~y.lag.1,subset=5:n)


#| echo: true
summary(model.lm)


#| echo: true
## Examine the first few fitted values from the nonparametric and
## linear AR(1) models
cbind(fitted(model.glp),fitted(model.lm))[1:10,]


#| echo: true
#| eval: false
## ## Fan & Yao's (1998) estimator, trim negative values of they occur
## ## Chen, Cheng & Peng's (2009) estimator is always positive but my
## ## Monte Carlo simulations show it is less efficient for this DGP
## library(np)
## set.seed(42)
## n <- 1000
## x <- sort(runif(n))
## sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
## dgp <- sin(2*pi*x)
## y <- dgp + rnorm(n,sd=sqrt(sigma.var))
## model <- npreg(y~x,regtype="ll")
## r <- residuals(model)**2
## ## Fan and Yao's (1998) estimator with trimming if needed
## var.fy <- fitted(npreg(r~x,regtype="ll"))
## var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## ## Chen, Cheng and Peng's (2009) estimator
## log.r <- log(r+1/n) ## Avoids log(0)
## V.hat <- fitted(npreg(log.r~x,regtype="ll"))
## d.hat <- 1/mean(r*exp(-V.hat))
## var.ccp <- exp(V.hat)/d.hat
## par(mfrow=c(1,2))
## plot(x,y,cex=.25,col="grey")
## lines(x,dgp,col=1,lty=1,lwd=2)
## lines(x,fitted(model),col=2,lty=2,lwd=2)
## legend("topleft",c("g(x)=sin(2 pi x)",
##                    "Nonparametric Estimate"),
##        col=1:2,
##        lty=1:2,
##        lwd=c(2,2),
##        bty="n")
## ylim=c(min(r),quantile(r,0.95))
## plot(x,r,ylim=ylim,cex=.25,col="grey")
## lines(x,sigma.var,col=1,lty=1,lwd=2)
## lines(x,var.fy,col=2,lty=2,lwd=2)
## lines(x,var.ccp,col=3,lty=4,lwd=2)
## legend("topleft",c("volatility=(2 pi x)^2/10",
##                    "Fan and Yao (1998)",
##                    "Chen et al. (2009)"),
##        col=1:3,
##        lty=c(1,2,4),
##        lwd=c(2,2,2),
##        bty="n")


#| label: fig-fanyaosim
#| fig-cap: Fan and Yao's (1998) and Chen et al.'s (2009) conditional variance estimators
#| echo: false
par(mfrow=c(1,2))
## Fan & Yao's (1998) estimator, trim negative values of they occur
## Chen, Cheng & Peng's (2009) estimator is always positive but my
## Monte Carlo simulations show it is less efficient for this DGP
set.seed(42)
n <- 1000
x <- sort(runif(n))
sigma.var <- 0.1*(.Machine$double.eps+2*pi*x)**2
dgp <- sin(2*pi*x)
y <- dgp + rnorm(n,sd=sqrt(sigma.var))
model <- npreg(y~x,regtype="ll")
r <- residuals(model)**2
## Fan and Yao's (1998) estimator with trimming if needed
var.fy <- fitted(npreg(r~x,regtype="ll"))
var.fy <- ifelse(var.fy<=0,.Machine$double.eps,var.fy)
## Chen, Cheng and Peng's (2009) estimator
log.r <- log(r+1/n) ## Avoids log(0)
V.hat <- fitted(npreg(log.r~x,regtype="ll"))
d.hat <- 1/mean(r*exp(-V.hat))
var.ccp <- exp(V.hat)/d.hat
plot(x,y,cex=.25,col="grey")
lines(x,dgp,col=1,lty=1,lwd=2)
lines(x,fitted(model),col=2,lty=2,lwd=2)
legend("topleft",c("g(x)=sin(2 pi x)", 
                   "Nonparametric Estimate"),
       col=1:2,
       lty=1:2,
       lwd=c(2,2),
       bty="n")
ylim=c(min(r),quantile(r,0.95))
plot(x,r,ylim=ylim,cex=.25,col="grey")
lines(x,sigma.var,col=1,lty=1,lwd=2)
lines(x,var.fy,col=2,lty=2,lwd=2)
lines(x,var.ccp,col=3,lty=4,lwd=2)
legend("topleft",c("volatility=(2 pi x)^2/10",
                   "Fan and Yao (1998)",
                   "Chen et al. (2009)"),
       col=1:3,
       lty=c(1,2,4),
       lwd=c(2,2,2),
       bty="n")


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## ## Compute the pooled first-step estimator
## model.pooled <- npreg(y~x,ckerorder=4)
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~x,ckerorder=2)


#| echo: true
#| eval: false
## ## Below is pseudo-code for generic data y, x (univariate) and
## ## panel id (univariate) - we shall modify and use this for real
## ## data but it demonstrates all the nuances of the approach
## library(np)
## model.fe <- npreg(y~x+factor(id))


#| echo: true
#| eval: false
## library(Ecdat)
## data(Airline)
## ?Airline
## library(np)
## options(np.messages=FALSE)
## attach(Airline)
## 
## ## Dependent variable is (log) cost, predictors (log) output, fuel price, load factor, year
## lcost <- as.numeric(log(cost))
## loutput <- as.numeric(log(output))
## lpf <- as.numeric(log(pf))
## lf <- as.numeric(lf)
## year <- ordered(year)
## airline <- factor(airline)
## 
## ## Airline specific fixed effects
## 
## model.fe <- npreg(lcost~loutput+lpf+lf+year+airline,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine",
##                   gradients=TRUE)
## 
## summary(model.fe)
## summary(model.fe)
## ## Plot partial means
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.fe,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)
## 
## ## Random effects
## 
## ## Compute the pooled first-step estimator
## 
## id <- airline
## n <- length(unique(airline))
## T <- length(unique(year))
## model.pooled <- npreg(lcost~loutput+lpf+lf+year,
##                       ckerorder=4,
##                       regtype="ll",
##                       bwmethod="cv.aic",
##                       ukertype="liracine",
##                       okertype="liracine")
## ## Compute Sigma.inv
## epsilon <- residuals(model.pooled)
## epsiloni.bar <- numeric()
## i <- sigmasq.v <- 0
## for(id in unique(id)) {
##   i <- i+1
##   epsiloni.bar[i] <- mean(epsilon[id==id])
##   sigmasq.v <- sigmasq.v + sum((epsilon[id==id]-epsiloni.bar[i])**2)
## }
## sigmasq.one <- (T/n)*sum(epsiloni.bar**2)
## sigmasq.v <- sigmasq.v/(n*(T-1))
## sigmasq.u <- (sigmasq.one-sigmasq.v)/T
## Sigma.inv <- (1/sigmasq.v)*diag(T) -
##   (sigmasq.u/sigmasq.v)/(sigmasq.v+sigmasq.u*T)*matrix(1,T,T)
## ## Compute Omega.inv.sqrt via Cholesky decomposition
## Omega.inv.sqrt <- chol(kronecker(diag(n),Sigma.inv))
## ## Compute c (for computing tau)
## c <- sigmasq.u/(sigmasq.u+sigmasq.v)
## ## Compute y-star
## y.star <- fitted(model.pooled) +
##   (1.5-c)*sqrt(sigmasq.u+sigmasq.v)*Omega.inv.sqrt%*%residuals(model.pooled)
## ## Compute the second-step estimator
## model.re <- npreg(y.star~loutput+lpf+lf+year,
##                   ckerorder=2,
##                   regtype="ll",
##                   bwmethod="cv.aic",
##                   ukertype="liracine",
##                   okertype="liracine")
## 
## summary(model.re)
## ## Plot partial means
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      common.scale=FALSE)
## ## Plot partial gradients
## plot(model.re,
##      plot.errors.method="bootstrap",
##      plot.errors.boot.num=25,
##      gradients=TRUE,
##      common.scale=FALSE)

#| echo: true
#| eval: true
library(np)
data(wage1)
model.lm <- lm(lwage~female+married+educ+tenure+exper+I(exper^2),data=wage1)
knitr::kable(coef(model.lm),col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.pl <- npplreg(lwage~female+married+educ+tenure|exper,
                    data=wage1)
knitr::kable(coef(model.pl),digits=7,col.names="Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.scoef <- npscoef(lwage~educ+tenure+exper+expersq|female+married,
                       data=wage1,
                       betas=TRUE)
summary(model.scoef)


#| echo: true
#| eval: true
## The smooth coefficients are vectors, one for each predictor in X (they "vary"
## with Z), so we compute the means of the columns of these coefficients, one for
## each X predictor
knitr::kable(colMeans(coef(model.scoef)),digits=7,col.names="Avg. Coefficient")


#| echo: true
#| eval: true
library(np)
data(wage1)
model.index <- npindex(lwage~educ+
                       tenure+
                       exper+
                       expersq+
                       female+
                       married,
                       method="ichimura",
                       data=wage1)
summary(model.index)


#| echo: true
#| eval: true
library(MASS)
data("birthwt")
model.logit <- glm(low~factor(smoke)+
                   factor(race)+
                   factor(ht)+
                   factor(ui)+
                   ordered(ftv)+
                   age+
                   lwt,
                   family=binomial(link=logit),data=birthwt)
cm.logit <- with(birthwt,table(low, ifelse(fitted(model.logit)>0.5, 1, 0)))
ccr.logit <- sum(diag(cm.logit))/sum(cm.logit)
knitr::kable(cm.logit)


#| echo: true
#| eval: true
library(MASS)
data(birthwt)
library(np)
model.index <- npindex(low~factor(smoke)+
                       factor(race)+
                       factor(ht)+
                       factor(ui)+
                       ordered(ftv)+
                       age+
                       lwt,
                       method="kleinspady",
                       data=birthwt)
cm.index <- with(birthwt,table(low, ifelse(fitted(model.index)>0.5, 1, 0)))
ccr.index <- sum(diag(cm.index))/sum(cm.index)
knitr::kable(cm.index)


require(np)
require(MASS)
## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(np)
options(np.messages=FALSE)
set.seed(42)
## This function takes a confusion matrix and formats it correctly if
## it is unbalanced and returns the CCR as well.
CM <- function(cm) {
  factor.values.eval <- colnames(cm)
  CM <- matrix(0,nrow(cm),nrow(cm))
  rownames(CM) <- rownames(cm)
  colnames(CM) <- rownames(cm)
  for(i in 1:ncol(cm)) CM[,(1:nrow(cm))[rownames(cm)==factor.values.eval[i]]] <- cm[,i]
  return(list(CM=CM,CCR=sum(diag(CM))/sum(CM)))
}
## Load the birthwt data
library(MASS)
data(birthwt)
## Create a data frame that has up to 4th order polynomials. We will
## use the BIC-optimal parametric model that will allow for
## interactions and some nonlinearity.
bwt <- with(birthwt,data.frame(low=factor(low),
                               race=factor(race),
                               smoke=factor(smoke),
                               ht=factor(ht),
                               ui=factor(ui),
                               ftv=ordered(ftv),
                               age=age,
                               agesq=age**2,
                               agecu=age**3,
                               agequ=age**4,
                               lwt=lwt,
                               lwtsq=lwt**2,
                               lwtcu=lwt**3,
                               lwtqu=lwt**4))
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), number of multistarts for bandwidth selection (nmulti)
## and number of train/eval splits (M)
n.eval <- 10
n <- nrow(bwt)
n.train <- n-n.eval
M <- 1000
## Create storage vectors
ccr.linear <- numeric(M)
ccr.linint <- numeric(M)
ccr.BIC <- numeric(M)
ccr.kernel <- numeric(M)
ccr.index <- numeric(M)
## Copy the full sample into the object train
train <- bwt
## Fit the parametric Logit model for the full sample.
##
## Linear 
logit.linear <- glm(low~
                    smoke+
                    race+
                    ht+
                    ui+
                    age+
                    lwt+
                    ftv,
                    family=binomial(link=logit),
                    data=train)
## Linear Logit model with interactions
logit.linint <- glm(low~
                    (smoke+
                     race+
                     ht+
                     ui+
                     age+
                     lwt+
                     ftv)^2,
                    family=binomial(link=logit),
                    data=train)
## BIC-optimal Logit model
logit.BIC <- glm(low ~ .,
                 family = binomial(link=logit),
                 data = train)
logit.BIC <- stepAIC(logit.BIC, ~ .^3,
                     trace=FALSE,
                     k=log(nrow(birthwt)))
## Klein spady estimator
model.index <- npindex(low~smoke+
                       race+
                       ht+
                       ui+
                       ftv+
                       age+
                       lwt,
                       method="kleinspady",
                       data=train)
## Get the bandwidths for the nonparametric model for the full sample.
bw <- npcdensbw(low~
                smoke+
                race+
                ht+
                ui+
                age+
                lwt+
                ftv,
                data=train)
## Apparent (in-sample) performance
ccr.app.linear <- with(train,CM(table(low,ifelse(predict(update(logit.linear),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.linint <- with(train,CM(table(low,ifelse(predict(update(logit.linint),
                                                 type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.BIC <- with(train,CM(table(low,ifelse(predict(update(logit.BIC),
                                              type="response",newdata=train)>0.5,1,0))))$CCR
ccr.app.index <- with(train,CM(table(low,ifelse(fitted(model.index)>0.5, 1, 0))))$CCR
ccr.app.kernel <- npconmode(bws=bw,newdata=train)$CCR.overall
## Conduct the M train/eval splits
for(m in 1:M) {
    ## Shuffle the data into independent training and evaluation samples.
    ii <- sample(1:n,replace=FALSE)
    train <- bwt[ii[1:n.train],]
    ## glm() can't deal with < all ftv cases
    while(length(unique(train$ftv))<length(unique(bwt$ftv))) {
        ii <- sample(1:n,replace=FALSE)
        train <- bwt[ii[1:n.train],]
    }
    eval <- bwt[ii[(n.train+1):n],]
    ## Extract the correct classification ratios for the independent
    ## evaluation data where we know the outcomes (update() refits the
    ## Logit model on train, the nonparametric model will
    ## automatically update taking train from the environment when it
    ## is called by predict()).
    ccr.linear[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linear),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.linint[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.linint),
                                                   type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.BIC[m] <- with(eval,CM(table(low,ifelse(predict(update(logit.BIC),
                                                type="response",newdata=eval)>0.5,1,0))))$CCR
    ccr.index[m] <- with(eval,CM(table(low,ifelse(predict(model.index,newdata=eval)>0.5,1,0))))$CCR
    ccr.kernel[m] <- npconmode(bws=bw,newdata=eval)$CCR.overall
}
## Conduct a paired t-test that the mean expected true CCR for each model 
## is equal versus the alternative that the kernel-based model has a significantly larger
## expected true CCR than the BIC-optimal Logit model.
p.linear <- t.test(ccr.kernel,ccr.linear,alternative="greater",paired=TRUE)$p.value
p.linint <- t.test(ccr.kernel,ccr.linint,alternative="greater",paired=TRUE)$p.value
p.BIC <- t.test(ccr.kernel,ccr.BIC,alternative="greater",paired=TRUE)$p.value
p.index <- t.test(ccr.kernel,ccr.index,alternative="greater",paired=TRUE)$p.value
p <- c(p.linear,p.linint,p.BIC,p.index,NA)
## Apparent performance
apparent <- c(ccr.app.linear,ccr.app.linint,ccr.app.BIC,ccr.app.index,ccr.app.kernel)
## Expected true performance
true <- c(mean(ccr.linear),mean(ccr.linint),mean(ccr.BIC),mean(ccr.index),mean(ccr.kernel))



foo <- data.frame(apparent,true,rank(-true),p)
rownames(foo) <- c("Par-Linear","Par-Linear-Int","Par-BIC","Semipar-KS","Nonpar")
colnames(foo) <- c("Apparent","Expected","Rank","$P$-value")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (higher values are preferred) and $P$-values from a test for equality of expected true performance. The cases considered are the kernel versus linear models, kernel versus linear with interaction models, kernel versus BIC-optimal models, and kernel versus the Klein and Spady models. Rejection of the null implies that the kernel-based model has significantly higher mean CCR on independent data.",
             digits=3,
             escape=FALSE)



## This code chunk illustrates the RP test (Revealed Performance)
## detailed in Racine and Parmeter (2012)
require(crs)
require(np)
data(wage1)
options(crs.messages=FALSE,np.messages=FALSE)
set.seed(42)
## Set the size of the evaluation data (n.eval) and the training data
## (n.train), and number of train/eval splits (M).
n.eval <- 10
n <- nrow(wage1)
n.train <- n-n.eval
M <- 1000
## Create storage vectors.
aspe.linear <- numeric(M)
aspe.quadint <- numeric(M)
aspe.mw <- numeric(M)
aspe.vc <- numeric(M)
aspe.pl <- numeric(M)
aspe.si <- numeric(M)
aspe.lc <- numeric(M)
aspe.ll <- numeric(M)
aspe.glp <- numeric(M)
## Copy the full sample into the object train.
train <- wage1
## Fit the parametric, semiparametric, and nonparametric models for 
## the full sample.
##
## Linear 
lm.linear <- lm(lwage~female+
                married+
                educ+
                exper+
                tenure,
                data=train)
## Linear with interactions and quadratic in experience (common
## specification).
lm.quadint <- lm(lwage~(female+
                        married+
                        educ+
                        exper+
                        I(exper^2)+
                        tenure)^2,
                data=train)
## Murphy-Welch quartic specification                
lm.mw <- lm(lwage~female+
            married+
            educ+
            exper+
            I(exper^2)+
            I(exper^3)+
            I(exper^4)+
            tenure,
            data=train)
## Varying coefficient
model.vc <- npscoef(lwage~educ+
                    tenure+
                    exper+
                    expersq|female+married,
                    data=train)
## Partially linear
model.pl <- npplreg(lwage~female+
                    married+
                    educ+
                    tenure|exper,
                    data=train)
## Single-index (Ichimura)
model.si <- npindex(lwage~female+
                    married+
                    educ+
                    tenure+
                    exper+
                    expersq,
                    method="ichimura",
                    data=train)
## Local constant
bw.lc <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="lc",
                 bwmethod="cv.aic",
                 data=train)
model.lc <- npreg(bws=bw.lc)
## Local linear
bw.ll <- npregbw(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 regtype="ll",
                 bwmethod="cv.aic",
                 data=train)
model.ll <- npreg(bws=bw.ll)
## Generalized local polynomial
ghat <- npglpreg(lwage~female+
                 married+
                 educ+
                 exper+
                 tenure,
                 cv.func="cv.aic",
                 data=train)
model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
## Apparent (in-sample) performance
aspe.app.linear <- with(train,mean((lwage-predict(lm.linear,newdata=train))^2))
aspe.app.quadint <- with(train,mean((lwage-predict(lm.quadint,newdata=train))^2))
aspe.app.mw <- with(train,mean((lwage-predict(lm.mw,newdata=train))^2))
aspe.app.vc <- with(train,mean((lwage-predict(model.vc,newdata=train))^2))
aspe.app.pl <- with(train,mean((lwage-predict(model.pl,newdata=train))^2))
aspe.app.si <- with(train,mean((lwage-predict(model.si,newdata=train))^2))
aspe.app.lc <- with(train,mean((lwage-predict(model.lc,newdata=train))^2))
aspe.app.ll <- with(train,mean((lwage-predict(model.ll,newdata=train))^2))
aspe.app.glp <- with(train,mean((lwage-predict(model.glp,newdata=train))^2))
## Conduct the M train/eval splits
for(m in 1:M) {
    ## We set the seed here to guarantee that the shuffles generated
    ## here and those in the Practitioner's Corner in Chapter 6 are
    ## identical.
    set.seed(m)
    ## Shuffle the data into independent training and evaluation
    ## samples.
    ii <- sample(1:n,replace=FALSE)
    train <- wage1[ii[1:n.train],]
    eval <- wage1[ii[(n.train+1):n],]
    ## Extract the APSEs for the independent evaluation data where we
    ## know the outcomes (update() refits the lm model on train).
    aspe.linear[m] <- with(eval,mean((lwage-predict(update(lm.linear),newdata=eval))^2))
    aspe.quadint[m] <- with(eval,mean((lwage-predict(update(lm.quadint),newdata=eval))^2))
    aspe.mw[m] <- with(eval,mean((lwage-predict(update(lm.mw),newdata=eval))^2))
    ## Calling the semi- and nonparametric functions with the existing
    ## bandwidth object re-estimates the model on the updated training data.
    model.vc.boot <- npscoef(bws=model.vc$bws)
    aspe.vc[m] <- with(eval,mean((lwage-predict(model.vc.boot,newdata=eval))^2))
    model.pl.boot <- npplreg(bws=model.pl$bw)
    aspe.pl[m] <- with(eval,mean((lwage-predict(model.pl.boot,newdata=eval))^2))
    model.si.boot <- npindex(bws=model.si$bws)
    aspe.si[m] <- with(eval,mean((lwage-predict(model.si.boot,newdata=eval))^2))
    model.lc <- npreg(bws=bw.lc)        
    aspe.lc[m] <- with(eval,mean((lwage-predict(model.lc,newdata=eval))^2))
    model.ll <- npreg(bws=bw.ll)
    aspe.ll[m] <- with(eval,mean((lwage-predict(model.ll,newdata=eval))^2))
    model.glp <- npglpreg(ghat$formula,cv="none",degree=ghat$degree,bws=ghat$bws,data=train)
    aspe.glp[m] <- with(eval,mean((lwage-predict(model.glp,newdata=eval))^2))
}
## Conduct a paired t-test that the mean expected true ASPE for each
## model is equal versus the alternative that the kernel has a
## significantly lower expected true ASPE than the MW-optimal lm
## model.
## LC versus parametric and semiparametric.
p.linear <- t.test(aspe.lc,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.lc,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.lc,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.lc,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.lc,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.lc,aspe.si,alternative="less",paired=TRUE)$p.value
p.lc <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## LL versus parametric and semiparametric.
p.linear <- t.test(aspe.ll,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.ll,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.ll,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.ll,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.ll,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.ll,aspe.si,alternative="less",paired=TRUE)$p.value
p.ll <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## GLP versus parametric and semiparametric.
p.linear <- t.test(aspe.glp,aspe.linear,alternative="less",paired=TRUE)$p.value
p.quadint <- t.test(aspe.glp,aspe.quadint,alternative="less",paired=TRUE)$p.value
p.mw <- t.test(aspe.glp,aspe.mw,alternative="less",paired=TRUE)$p.value
p.vc <- t.test(aspe.glp,aspe.vc,alternative="less",paired=TRUE)$p.value
p.pl <- t.test(aspe.glp,aspe.pl,alternative="less",paired=TRUE)$p.value
p.si <- t.test(aspe.glp,aspe.si,alternative="less",paired=TRUE)$p.value
p.glp <- c(p.linear,p.quadint,p.mw,p.vc,p.pl,p.si,NA,NA,NA)
## Apparent performance.
apparent <- c(aspe.app.linear,aspe.app.quadint,aspe.app.mw,
              aspe.app.vc,aspe.app.pl,aspe.app.si,
              aspe.app.glp,aspe.app.ll,aspe.app.lc)
## Expected true performance.
true <- c(mean(aspe.linear),mean(aspe.quadint),mean(aspe.mw),
          mean(aspe.vc),mean(aspe.pl),mean(aspe.si),
          mean(aspe.glp),mean(aspe.ll),mean(aspe.lc))



foo <- data.frame(apparent,true,rank(true),p.glp,p.ll,p.lc)
rownames(foo) <- c("Par-Linear","Par-Quad-Int","Par-Murphy-Welch",
                   "Semipar-VC","Semipar-PL","Semipar-SI",
                   "Nonpar-GLP","Nonpar-LL","Nonpar-LC")
colnames(foo) <- c("Apparent","Expected","Rank",
                   "$P$-GLP","$P$-LL","$P$-LC")
knitr::kable(foo,
             caption="Apparent versus expected true model performance (lower values are preferred) and $P$-values from a test for equality of expected true performance for the kernel versus the parametric and semiparametric models (rejection of the null implies the kernel model has significantly lower mean ASPE on independent data).",
             digits=4,
             escape=FALSE)


