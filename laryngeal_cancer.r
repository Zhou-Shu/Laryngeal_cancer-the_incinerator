library(spatstat)
cases <- split(chorley)$larynx
controls <- split(chorley)$lung
incin <- chorley.extra$incin
kernel_lung <- density(controls,sigma =bw.diggle(chorley), eps=0.1, positive=TRUE)
kernel_lung <- eval.im(pmax(kernel_lung, 1e-10))
dfun <- function(x,y) { sqrt((x-incin$x)^2 + (y-incin$y)^2) }
dinc <- as.im(dfun, chorley$window)
raisin <- function(x,y, exp_theta1, theta2) {1+exp_theta1 * exp( - theta2 * dfun(x,y))}

fit1 <- ippm(cases ~ offset(log(kernel_lung) + log(raisin)),start=list(exp_theta1=5, theta2=1))
fit2 <- update(chorleyDfit, . ~offset(log(kernel_lung)))

plot(kernel_lung,main="Lung Kernel")
plot(dinc,main="Distance to incinerate")
plot(fit1,main="Alternative model")
plot(fit2,main="Null model")
anova(fit2, fit1, test="LRT")
