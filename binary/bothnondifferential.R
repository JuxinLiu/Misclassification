# both Y and X subject to nondifferential misclassification errors

# true values of P(X=1), beta0, beta1, SN, SP, and D parameters

p1.true <- 0.4
beta0.true <- 1
beta1.true <- 1
yp <- c(1, 1, 0, 0)
xp <- c(1, 0, 1, 0)
# P(Y=y,X=x) in order of 11, 10, 01, 00
pi.true <- exp(yp * (beta0.true + beta1.true * xp))/(1 + exp(beta0.true + beta1.true * xp)) * p1.true^xp * (1 - p1.true)^(1 - xp)

# function that returns D values that produce max and min \delta
DD <- function(para) {
  eps = 1e-04
  SNY <- para[1]
  SNX <- para[2]
  SPY <- para[3]
  SPX <- para[4]
  L1 <- max(-SNY * SNX, -(1 - SNY) * (1 - SNX))
  U1 <- min((1 - SNY) * SNX, SNY * (1 - SNX))
  L2 <- max(-SNY * SPX, -(1 - SNY) * (1 - SPX))
  U2 <- min((1 - SNY) * SPX, SNY * (1 - SPX))
  L3 <- max(-SPY * SNX, -(1 - SPY) * (1 - SNX))
  U3 <- min((1 - SPY) * SNX, SPY * (1 - SNX))
  L4 <- max(-SPY * SPX, -(1 - SPY) * (1 - SPX))
  U4 <- min((1 - SPY) * SPX, SPY * (1 - SPX))
  obj <- function(D) {
    ans <- abs(sum(pi.true * D * c(1, -1, -1, 1)))
  }
  obj_minus <- function(D) {
    ans <- -abs(sum(pi.true * D * c(1, -1, -1, 1)))
  }
  Dlow <- optim(c(0, 1e-05, 1e-05, 1), obj, lower = c(L1, L2, L3, L4) + eps, upper = c(U1, U2, U3, U4) - eps, method = "L-BFGS-B")
  Dhigh <- optim(c(0, 1e-05, 1e-05, 1), obj_minus, lower = c(L1, L2, L3, L4) + eps, upper = c(U1, U2, U3, U4) - eps, method = "L-BFGS-B")
  ans <- list(D_l = round(c(Dlow$par, Dlow$val), 3), D_h = round(c(Dhigh$par, -Dhigh$val), 3))
  print(ans)
}

#### the columns in order of SNY, SNX, SPY, SPX
SNSP <- matrix(c(0.9, 0.9, 0.9, 0.9, 0.8, 0.8, 0.8, 0.8, 0.7, 0.7, 0.7, 0.7), nrow = 3, byrow = TRUE)

allcases <- matrix(NA, nrow = 6, ncol = 8)

for (s in 1:3) {
  DDtmp <- DD(SNSP[s, ])
  allcases[(s - 1) * 2 + 1, ] <- c(SNSP[s, ], DDtmp$D_l[1:4])
  allcases[2 * s, ] <- c(SNSP[s, ], DDtmp$D_h[1:4])
}

### data generation

para <- allcases[2, ]  # SN high,SP high,D_high

SNY.t <- SNX.t <- para[1]
SPY.t <- SPX.t <- para[3]
DD.t <- para[-(1:4)]

Ct <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  # i index observed values j index true values
  for (j in 1:4) {
    ystar <- 2 - ceiling(i/2)
    xstar <- ifelse(floor(i/2) == ceiling(i/2), 0, 1)
    y <- 2 - ceiling(j/2)
    x <- ifelse(floor(j/2) == ceiling(j/2), 0, 1)
    Ct[i, j] <- (SNY.t^ystar * (1 - SNY.t)^(1 - ystar))^y * ((1 - SPY.t)^ystar * SPY.t^(1 - ystar))^(1 - y) * (SNX.t^xstar * (1 -
                                                                                                                                SNX.t)^(1 - xstar))^x * ((1 - SPX.t)^xstar * SPX.t^(1 - xstar))^(1 - x) + DD.t[j] * (-1)^(ifelse(ystar == y, 1, 0) + ifelse(xstar ==
                                                                                                                                                                                                                                                              x, 1, 0))
    }  # end of j loop
}  # end of i loop

Ct <- ifelse(Ct < 1e-15, 0, Ct)

# validation data (ystar,xstar|y,x)(y,x) by columns (y,x) in order of 11, 10, 01, 00
pjoint.true <- as.vector(t(t(Ct) * pi.true))

Nv <- 1000
Nm <- 9000

set.seed(123)
validation <- as.vector(rmultinom(1, Nv, pjoint.true))
pstar.true <- Ct %*% pi.true[1:4]
main <- as.vector(rmultinom(1, Nm, pstar.true))

# validation data (ystar,xstar|y,x)(y,x) by columns (y,x) in order of 11, 10, 01, 00

xlabel <- ylabel <- 1:0
all <- as.data.frame(expand.grid(xlabel, ylabel, xlabel, ylabel))
names(all) <- c("xstar", "ystar", "x", "y")

all$count <- validation
SNY_v <- sum(all$count[all$ystar == 1 & all$y == 1])/sum(all$count[all$y == 1])
SPY_v <- sum(all$count[all$ystar == 0 & all$y == 0])/sum(all$count[all$y == 0])
SNX_v <- sum(all$count[all$xstar == 1 & all$x == 1])/sum(all$count[all$x == 1])
SPX_v <- sum(all$count[all$xstar == 0 & all$x == 0])/sum(all$count[all$x == 0])

D0 <- rep(NA, 1)
D0[1] <- all$count[all$ystar == 1 & all$xstar == 1 & all$y == 1 & all$x == 1]/sum(all$count[all$y == 1 & all$x == 1]) - SNY_v * SNX_v
D0[2] <- all$count[all$ystar == 1 & all$xstar == 0 & all$y == 1 & all$x == 0]/sum(all$count[all$y == 1 & all$x == 0]) - SNY_v * SPX_v
D0[3] <- all$count[all$ystar == 0 & all$xstar == 1 & all$y == 0 & all$x == 1]/sum(all$count[all$y == 0 & all$x == 1]) - SPY_v * SNX_v
D0[4] <- all$count[all$ystar == 0 & all$xstar == 0 & all$y == 0 & all$x == 0]/sum(all$count[all$y == 0 & all$x == 0]) - SPY_v * SPX_v

px0 <- sum(all$count[all$x == 1])/sum(all$count)

logit <- function(t) {
  ans <- log(t) - log(1 - t)
}
beta0_ini <- logit(sum(all$count[all$y == 1 & all$x == 0])/sum(all$count[all$x == 0]))
beta1_ini <- logit(sum(all$count[all$y == 1 & all$x == 1])/sum(all$count[all$x == 1])) - beta0_ini

library(coda)
library(runjags)

### data

a1 <- 2.7
b1 <- 9  # 95% ET-CI (0.05,0.5)

data <- list(validation = validation, main = main, Nm = Nm, Nv = Nv, a1 = a1, b1 = b1)


### indepdent model

indmodel <- "
model{
    for (i in 1:4) {
    # i index observed values j index true values
    for (j in 1:4) {
        ystar[i, j] <- 2 - round(i/2)
        xstar[i, j] <- ifelse(round(i/2) == trunc(i/2), 0, 1)
        y[i, j] <- 2 - round(j/2)
        x[i, j] <- ifelse(round(j/2) == trunc(j/2), 0, 1)
        # C[i,j]=P(Ystar=ystar,Xstar=xstar|Y=y,X=x); i indexing the pair (ystar,xstar) and j indexing the pair (y,x)
        C[i, j] <- (SNY^ystar[i, j] * (1 - SNY)^(1 - ystar[i, j]))^y[i, j] * ((1 - SPY)^ystar[i, j] * SPY^(1 - ystar[i, j]))^(1 - y[i,
            j]) * (SNX^xstar[i, j] * (1 - SNX)^(1 - xstar[i, j]))^x[i, j] * ((1 - SPX)^xstar[i, j] * SPX^(1 - xstar[i, j]))^(1 - x[i,
            j])
    }  # end of j loop
}  # end of i loop

# validation data
validation[1:16] ~ dmulti(pjoint, Nv)
p[1] <- exp(beta0 + beta1)/(1 + exp(beta0 + beta1)) * p1  #[1,1]
p[2] <- exp(beta0)/(1 + exp(beta0)) * (1 - p1)  #[1,0]
p[3] <- 1/(1 + exp(beta0 + beta1)) * p1  #[0,1]
p[4] <- 1/(1 + exp(beta0)) * (1 - p1)  #[0,0]

# p[1:4]<-exp(yp*(beta0+beta1*xp))/(1+exp(beta0+beta1*xp))*exp(xp*log(p1)+(1-xp)*log(1-p1))
for (s in 1:4) {
    pjoint[(4 * (s - 1) + 1):(4 * s)] <- C[, s] * p[s]
}

# main data
main[1:4] ~ dmulti(pstar[1:4], Nm)
pstar[1:4] <- C[1:4, 1:4] %*% p[1:4]

# priors
p1 ~ dbeta(a1, b1)
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)

SNX ~ dunif(1 - SPX, 1)
SPX ~ dunif(0.5, 1)
SNY ~ dunif(1 - SPY, 1)
SPY ~ dunif(0.5, 1)
}
"

monitor <- c("SNY", "SNX", "SPY", "SPX", "beta0", "beta1")  # names of all the parameters
# note both p and cov are in the ordering of 11,10,01,00

inits1 <- list(SNY = SNY_v, SPY = SPY_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, beta0 = -1, beta1 = 2, .RNG.name = "base::Super-Duper", .RNG.seed = 1)
inits2 <- list(SNY = SNY_v, SPY = SPY_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, beta0 = 0, beta1 = 1, .RNG.name = "base::Wichmann-Hill",
               .RNG.seed = 2)

opt <- run.jags(model = indmodel, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000, method = "rjags",
                inits = list(inits1, inits2))

out <- summary(opt)[, c(2, 4, 5, 11)]

########### naive model

naive <- "
model{
    for (i in 1:4) {
    # i index observed values j index true values
    for (j in 1:4) {
        ystar[i, j] <- 2 - round(i/2)
        xstar[i, j] <- ifelse(round(i/2) == trunc(i/2), 0, 1)
        y[i, j] <- 2 - round(j/2)
        x[i, j] <- ifelse(round(j/2) == trunc(j/2), 0, 1)
        # C[i,j]=P(Ystar=ystar,Xstar=xstar|Y=y,X=x); i indexing the pair (ystar,xstar) and j indexing the pair (y,x)
        C[i, j] <- exp(ystar[i, j] * (beta0 + xstar[i, j] * beta1))/(1 + exp(beta0 + xstar[i, j] * beta1)) * p1^(xstar[i, j]) * (1 - p1)^(1 -
            xstar[i, j])
    }  # end of j loop
}  # end of i loop

# validation data
validation[1:16] ~ dmulti(pjoint[1:16], Nv)
p[1] <- exp(beta0 + beta1)/(1 + exp(beta0 + beta1)) * p1  #[1,1]
p[2] <- exp(beta0)/(1 + exp(beta0)) * (1 - p1)  #[1,0]
p[3] <- 1/(1 + exp(beta0 + beta1)) * p1  #[0,1]
p[4] <- 1/(1 + exp(beta0)) * (1 - p1)  #[0,0]
for (s in 1:4) {
    pjoint[(4 * (s - 1) + 1):(4 * s)] <- C[, s] * p[s]
}

# main data
main[1:4] ~ dmulti(p[1:4], Nm)  # no misclassification in main data

# priors
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
p1 ~ dbeta(a1, b1)
}
"
data_naive <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, a1 = a1, b1 = b1)
monitor <- c("beta0", "beta1")  # names of all the parameters
# note both p and cov are in the ordering of 11,10,01,00

inits1 <- list(beta0 = -1, beta1 = 2, p1 = 0.2, .RNG.name = "base::Super-Duper", .RNG.seed = 1)
inits2 <- list(beta0 = 0, beta1 = 1, p1 = 0.2, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2)

opt_naive <- run.jags(model = naive, monitor = monitor, data = data_naive, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000,
                      method = "rjags", inits = list(inits1, inits2))

out_naive <- summary(opt_naive)[, c(2, 4, 5, 11)]

################# rightmodel

data_cov <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, a1 = a1, b1 = b1)

covmodel <- "
model{
    for (i in 1:4) {
    # i index observed values j index true values
    for (j in 1:4) {
        ystar[i, j] <- 2 - round(i/2)
        xstar[i, j] <- ifelse(round(i/2) == trunc(i/2), 0, 1)
        y[i, j] <- 2 - round(j/2)
        x[i, j] <- ifelse(round(j/2) == trunc(j/2), 0, 1)
        # C[i,j]=P(Ystar=ystar,Xstar=xstar|Y=y,X=x); i indexing the pair (ystar,xstar) and j indexing the pair (y,x)
        C[i, j] <- (SNY^ystar[i, j] * (1 - SNY)^(1 - ystar[i, j]))^y[i, j] * ((1 - SPY)^ystar[i, j] * SPY^(1 - ystar[i, j]))^(1 - y[i,
            j]) * (SNX^xstar[i, j] * (1 - SNX)^(1 - xstar[i, j]))^x[i, j] * ((1 - SPX)^xstar[i, j] * SPX^(1 - xstar[i, j]))^(1 - x[i,
            j]) + D[j] * (-1)^(ifelse(ystar[i, j] == y[i, j], 1, 0) + ifelse(xstar[i, j] == x[i, j], 1, 0))
    }  # end of j loop
}  # end of i loop

# validation data
validation[1:16] ~ dmulti(pjoint[1:16], Nv)
p[1] <- exp(beta0 + beta1)/(1 + exp(beta0 + beta1)) * p1  #[1,1]
p[2] <- exp(beta0)/(1 + exp(beta0)) * (1 - p1)  #[1,0]
p[3] <- 1/(1 + exp(beta0 + beta1)) * p1  #[0,1]
p[4] <- 1/(1 + exp(beta0)) * (1 - p1)  #[0,0]
for (s in 1:4) {
    pjoint[(4 * (s - 1) + 1):(4 * s)] <- C[, s] * p[s]
}

# main data
main[1:4] ~ dmulti(pstar[1:4], Nm)
pstar[1:4] <- C[1:4, 1:4] %*% p[1:4]

# priors
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01)
p1 ~ dbeta(a1, b1)
SNX ~ dunif(1 - SPX, 1)
SPX ~ dunif(0.5, 1)
SNY ~ dunif(1 - SPY, 1)
SPY ~ dunif(0.5, 1)

L[1] <- max(-SNY * SNX, -(1 - SNY) * (1 - SNX))
U[1] <- min((1 - SNY) * SNX, SNY * (1 - SNX))
L[2] <- max(-SNY * SPX, -(1 - SNY) * (1 - SPX))
U[2] <- min((1 - SNY) * SPX, SNY * (1 - SPX))
L[3] <- max(-SPY * SNX, -(1 - SPY) * (1 - SNX))
U[3] <- min((1 - SPY) * SNX, SPY * (1 - SNX))
L[4] <- max(-SPY * SPX, -(1 - SPY) * (1 - SPX))
U[4] <- min((1 - SPY) * SPX, SPY * (1 - SPX))

# priors for D parameters
for (j in 1:4) {
    Dt[j] ~ dunif(0, 1)
    D[j] <- L[j] + (U[j] - L[j]) * Dt[j]
  }
}#end of model
"

monitor <- c("SNY", "SNX", "SPY", "SPX", "D", "beta0", "beta1")  # names of all the parameters
# note both p and cov are in the ordering of 11,10,01,00

inits1 <- list(SNY = SNY_v, SPY = SPY_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, Dt = rep(0.8, 4), beta0 = -1, beta1 = 2, .RNG.name = "base::Super-Duper",
               .RNG.seed = 1)
inits2 <- list(SNY = SNY_v, SPY = SPY_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, Dt = rep(0.5, 4), beta0 = 0, beta1 = 1, .RNG.name = "base::Wichmann-Hill",
               .RNG.seed = 2)

opt_cov <- run.jags(model = covmodel, monitor = monitor, data = data_cov, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000, method = "rjags",
                    inits = list(inits1, inits2))

out_cov <- summary(opt_cov)[, c(2, 4, 5, 11)]