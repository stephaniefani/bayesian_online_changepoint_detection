#Bayesian Online Change Point Detection
#Adams & Mackay (2007)
#Using Exponential Conjugate (Normal - Inverse Gamma Distribution)

#Generate dx
dx <- rep(1,100)
for (i in 1:100) {
  if (i<= 40) {
    dx[i] <- rnorm(1,500,5) #normal distribution with mean 500
  }
  else {
    dx[i] <- rnorm(1,1200,5) #normal distribution with mean 1200
  }
}

#Plot x
plot(dx,
     main="Plot Harga Buka Saham ANTM",
     ylab="harga",
     type="l",
     col="blue")

#Let x ~ N(mean, var), with parameters: mean ~ N(u,var/w) and var ~ Inverse Gamma(a,b) 

#Prior hyperparameter 
#Byrd, Michael; Nghiem, Linh; Cao, Jing. 2017.Lagged Exact Bayesian Online Changepoint Detection with Parameter Estimation
u0 <- 500
w0 <- 10^-4
a0 <- 1
b0 <- 10^-5

#prior hyperparameter tuning is needed to improve the model

#loop
#Initiate r (run length)
r <- 0
max <- rep(0,100)
i <- 1

while (i <= 100) {
  if (r==0)  #if changepoint occurs
  {
    #Initiate prior hyperparameter
    u <- u0
    w <- w0
    a <- a0
    b <- b0
    
    #Observe data
    x <- dx[i]
    
    #Predictive Distribution P(xt |xt-1) ~  T (unstandardized student-t)
    k <- u
    l <- (b*(w+1))/a*w
    v <- 2*a #degree of freedom
    
    p <- rep(0,1)
    p[r+2] <- (gamma((v+1)/2) / (sqrt(pi*v*l)*gamma(v/2))) * ((1+1/v*((x-k)^2)/l)^(-(v+1)/2))
    
    #Hazard Function
    library(survival)
    library(survminer)
    set.seed(7)
    time<-seq(1, 20, by=1)
    stat <- sample(c(0,1), replace=TRUE, size=20)
    s1 <- Surv(time,stat)
    model1<-survreg(s1 ~ 1, dist="exponential")
    summary(model1)
    h <- exp(-(summary(model1)$table[,1]))
    h <- 1/40 #set lambda
    
    #h = 0.04285714
    
    #Growth / Change point Function
    g <- matrix(c(p[r+2], -p[r+2], 0, p[r+2]),byrow=TRUE,nrow=2,ncol=2) %*% matrix(c(1,h),byrow=TRUE,nrow=2,ncol=1)
    
    #Run length distribution
    d <- 1/sum(g) * g
    
    #Finding max 
    max[i] <- which.max(d)
    
    r <- r+1
  } else #if change point not occurs
  {
    #Posterior hyperparameter (updating using data/ Bayesian updating)
    u <- (x + w*u)/(w + 1)
    w <- w + 1
    a <- a + 1/2
    b <- b + (w*((x-u)^2))/(2*(w+1)) 
    
    #Observe data
    x <- dx[i]
    
    #Predictive Distribution P(xt |xt-1)
    k <- u
    l <- (b*(w+1))/a*w
    v <- 2*a #degree of freedom
    
    p[r+2] <- (gamma((v+1)/2) / (sqrt(pi*v*l)*gamma(v/2))) * ((1+1/v*((x-k)^2)/l)^(-(v+1)/2))
    
    #Growth / Change point Function
    
    g <- (diag(rev(p)*(1-h)) + matrix(c(rep(0,((r+1)*(r+2))),rev(p)*h),byrow=TRUE,nrow=r+2,ncol=r+2)) %*% matrix(c(g,0),byrow=TRUE,nrow=r+2,ncol=1)
    
    #Run length distribution
    d <- 1/sum(g) * g
    
    #Finding max 
    max[i] <- which.max(d)
    
    if ((r+3)-max[i] < 0.5*(r+2))
    {
      i  <- i - ((r+3)-max[i])
      r <- 0
      max[i-1] <- 0
    }else
    {
      r <- r+1
    }
  }
  print (r)
  print (i)
  i <- i+1
}

run <- rep(1,100)
for (j in 2:100)
{
  if (max[j] == 1)
  {
    run[j] <- run[j-1]+1
  } else
  {
    run[j] <- 0
  }
}

t <- seq(1,100,by=1)
plot(t,run,
     main="Plot Run Length",
     ylab="run length",
     type="l",
     col="black")


plot(t,dx,
     main="Plot Harga Buka Saham ANTM",
     type="l",
     col="blue",
     ylim = c(0,1200))
lines(t,run,col="black")


