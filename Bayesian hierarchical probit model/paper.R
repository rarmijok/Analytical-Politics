paper <- read.csv("paperdata.csv")
hist(paper$NSE,main="Distribution by social class")
plot(paper$xpos,paper$vote, main="Figure 1: GLM Probit model for probability of vote choice", xlab="Ideology",ylab="Probability of voting for Piñera")
abline(h = c(0,1), lty = 3)
# (paper$xpos,paper$Voto1)

probit.1 <- glm(paper$vote ~ paper$xpos, family = binomial(link = "probit"))
curve(pnorm(coef(probit.1)[1] + coef(probit.1)[2]*x),
      add = TRUE)
      
library(rjags)
library(coda)

probit.model <- jags.model("papermodel.bug",
                             data=list( "class"=paper$NSE, "y"=paper$vote, "x"=paper$xpos))

update(probit.model, n.iter=5000)
probit.samples <- coda.samples(model = probit.model,
                              variable.names = c("alpha", "beta"),
                              n.iter = 100000)


plot(probit.samples[[1]][,1:3])
plot(probit.samples[[1]][,4:6])


acf(probit.samples[[1]][,1], main="ACF for upper class intercept")
acf(probit.samples[[1]][,2], main="ACF for middle class intercept")
acf(probit.samples[[1]][,3], main="ACF for lower class intercept")
acf(probit.samples[[1]][,4], main="ACF for upper class slope")
acf(probit.samples[[1]][,5], main="ACF for middle class slope")
acf(probit.samples[[1]][,6], main="ACF for lower class slope")
effectiveSize(probit.samples)


geweke.diag(probit.samples, frac1=.1, frac2=.5)

plot(paper$xpos,paper$vote, main="Figure : Hierarchical probit model of vote choice", xlab="Ideology",ylab="Probability of voting for Piñera")
abline(h = c(0,1), lty = 3)
curve(pnorm(coef(probit.1)[1] + coef(probit.1)[2]*x), add = TRUE, lty=1)
curve(pnorm(mean(probit.samples[[1]][,1])+mean(probit.samples[[1]][,4])*x), add = TRUE, lty = 3)
curve(pnorm(mean(probit.samples[[1]][,2])+mean(probit.samples[[1]][,5])*x), add = TRUE, lty = 4)
curve(pnorm(mean(probit.samples[[1]][,4])+mean(probit.samples[[1]][,6])*x), add = TRUE, lty = 2)
legend(-2.2, 0.9, c("GLM","Upper class","Middle class","Lower class"), lty=c(1,3,4,2))


