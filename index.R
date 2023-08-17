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


