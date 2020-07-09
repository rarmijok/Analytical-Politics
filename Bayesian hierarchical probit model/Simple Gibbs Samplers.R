

#### MVN Gibbs Sampler Example



plot(sim1$x1, sim1$x2, type = "s")

plot(sim1$x1, sim1$x2, type = "s", xlim = c(-3,3), ylim = c(-3,3))

## trying three other samplers started at different values
sim2 <- mvn.gibbs(n = 1000, r = .7, inits = c(-10,-10), mu = c(0,0))
sim3 <- mvn.gibbs(n = 1000, r = .7, inits = c(10,-10), mu = c(0,0))
sim4 <- mvn.gibbs(n = 1000, r = .7, inits = c(-10,10), mu = c(0,0))

plot(sim1$x1, sim1$x2, xlim = c(-11,11), ylim = c(-11,11), type = "s")
lines(sim2$x1, sim2$x2, col = "red", type = "s")
lines(sim3$x1, sim3$x2, col = "green", type = "s")
lines(sim4$x1, sim4$x2, col = "blue", type = "s")


## looking at histograms
par(mfrow = c(1,2))
hist(sim1$x1)
hist(sim1$x2)
## what's wrong here?

## now, doing it omitting the first few draws
hist(sim1$x1[-(1:5)])
hist(sim1$x2[-(1:5)])



## looking at the samples of each parameter by itself
par(mfrow = c(2,1))
plot(sim1$x1, type = "l")
plot(sim1$x2, type = "l")


## now, looking at multiple chains on the same traceplot

par(mfrow = c(1,2))
plot(sim1$x1, ylim = c(-10, 10), type = "l")
lines(sim2$x1, col = "blue")
lines(sim3$x1, col = "red")
lines(sim4$x1, col = "green")

plot(sim1$x2, ylim = c(-10, 10), type = "l")
lines(sim2$x2, col = "blue")
lines(sim3$x2, col = "red")
lines(sim4$x2, col = "green")

 
## Calculating probability x1 > x2
mean(sim1$x1 > sim1$x2)


##################################################
## Gibbs sampler for 
## normal model with unknown mean and variance:
##########################


library(MCMCpack) # to get rinvgamma()

norm.samp <- function(m = 0, s2 = 1,
                      alpha = 1, beta = 1,
                      x,
                      startvals = c(0, 1), n.samples = 100){

  n <- length(x)
  xbar <- mean(x)

  # defining matrix for samples
  samples <- matrix(NA, nrow = n.samples, ncol = 2)
  colnames(samples) <- c("mu", "sigma2")

  # putting starting values as first row of samples matrix
  samples[1, ] <- startvals

  # looping for all samples
  for(i in 2:n.samples){
  
    # sampling from p(mu | x, sigma2)
    samples[i, 1] <- rnorm(n = 1,
                           mean = (m/s2 + n*xbar/samples[i-1, 2]) /
                           (1/s2 + n/samples[i-1, 2]),
                           sd = sqrt( 1 / (1/s2 + n/samples[i-1, 2])))
  
    xtilde <- 1/n * sum((x - samples[i-1, 1])^2)
    
    # sampling from p(sigma2 | x, mu)
    samples[i, 2] <- rinvgamma(n = 1,
                               alpha + n/2,
                               beta + n/2 * xtilde)
  }

  return(samples)
}



nsamp1 <- norm.samp(x = c(0,1,0,1,2,-1), n.samples = 10000)

plot(nsamp1[,1], nsamp1[,2], type = "s")

plot(density(nsamp1[,1]))
plot(density(nsamp1[,2]))

############################################################
############################################################
## Hierarchical Normal Model for average state ideologies ##
############################################################
############################################################

## generating fake dataset:

state.theta <- rnorm(10, mean = 50, sd = 5)
state.N <- c(3,5,10,13,22,36,55,66,142,154)

ideol <- c()
state <- c()
for(i in 1:10){
    foo1 <- rnorm(state.N[i], state.theta[i], 10)
    foo2 <- rep(i, state.N[i])
    ideol <- c(ideol, foo1)
    state <- c(state, foo2)
  }

#hnorm.data <- dget("Hierarchical Normal Model Data.txt")
#attach(hnorm.data)



## looking at raw state means
tapply(ideol, state, mean)
plot(tapply(ideol, state, mean) ~ state.N)
abline(h = mean(ideol))


theta.samp <- function(mu, sigmasq, tausq){
  ybar <- tapply(ideol, state, mean)
  theta <- rep(NA, 10)
  for(j in 1:10){
    mean <- ((1/tausq)*mu + (state.N[j]/sigmasq)*ybar[j])/
      ((1/tausq) + (state.N[j]/sigmasq))
    var <- 1 / ((1/tausq) + (state.N[j]/sigmasq))
    theta[j] <- rnorm(1, mean = mean, sd = sqrt(var))
  }
  return(theta)
}


mu.samp <- function(theta, tausq){
  thetabar <- mean(theta)
  mean <- ((1/20^2)*50 + (10/tausq)*thetabar) /
    ((1/20^2) + (10/tausq))
  var <- 1 / ((1/20^2) + (10/tausq))
  mu <- rnorm(1, mean = mean, sd = sqrt(var))
  return(mu)
}


library(MCMCpack)

sigmasq.samp <- function(theta){
  theta.state <- theta[state]
  y.twidle <- mean((ideol - theta.state)^2)
  sigmasq <- rinvgamma(1, .2 + length(state)/2, 5 + length(state)/2 * y.twidle)
  return(sigmasq)
}


tausq.samp <- function(theta, mu){
  theta.twidle <- mean((theta - mu)^2)
  tausq <- rinvgamma(1, .2 + 10/2, 5 + 10/2 * theta.twidle)
  return(tausq)
}


nsamples <- 1000
SAMPLES <- matrix(NA, nsamples, 13)
colnames(SAMPLES) <- c("theta1", "theta2", "theta3", "theta4", "theta5",
                       "theta6", "theta7", "theta8", "theta9", "theta10",
                       "mu", "sigmasq", "tausq")

SAMPLES[1, 11] <- rnorm(1, 50, sqrt(20))
SAMPLES[1, 12] <- runif(1, 0, 20)
SAMPLES[1, 13] <- runif(1, 0, 20)
SAMPLES[1, 1:10] <- rnorm(10, SAMPLES[1, 11], sqrt(SAMPLES[1, 13]))

for(i in 2:nsamples){
  SAMPLES[i, 1:10] <- theta.samp(mu = SAMPLES[i-1, 11],
                                 sigmasq = SAMPLES[i-1, 12],
                                 tausq = SAMPLES[i-1, 13])
  SAMPLES[i, 11] <- mu.samp(theta = SAMPLES[i, 1:10],
                            tausq = SAMPLES[i-1, 13])
  SAMPLES[i, 12] <- sigmasq.samp(theta = SAMPLES[i, 1:10])
  SAMPLES[i, 13] <- tausq.samp(theta = SAMPLES[i, 1:10],
                               mu = SAMPLES[i, 11])
}



par(mfrow = c(3, 1))
plot(SAMPLES[,11], type = "l", main = "Mu Traceplot")
plot(SAMPLES[,12], type = "l", main = "Sigma^2 Traceplot")
plot(SAMPLES[,13], type = "l", main = "Tau^2 Traceplot")

par(mfrow = c(5,2))
for(i in 1:10){
  plot(SAMPLES[, i], type = "l", main = paste("Traceplot Theta", i),
       ylim = c(30,65))
}

par(mfrow = c(5,2))
for(i in 1:10){
  plot(density(SAMPLES[, i]), main = paste("Density of Theta", i),
       xlim = c(30, 65))
}

## looking at state mean estimates alongside raw state means
tapply(ideol, state, mean)
plot(tapply(ideol, state, mean) ~ state.N)
points(apply(SAMPLES[ , 1:10], 2, mean) ~ state.N, col = "grey")
for(i in 1:10){
  segments(state.N[i], tapply(ideol, state, mean)[i],
           state.N[i], mean(SAMPLES[ , i]),
           col = "grey", cex = .3, lty = 3)
}
abline(h = mean(ideol))
