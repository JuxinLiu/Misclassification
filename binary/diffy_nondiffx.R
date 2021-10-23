# Y subject to differential errors and X subject to nondifferential errors

# true values of P(X=1), beta0, beta1, SN, SP, and D parameters
p1.true <- 0.4
beta0.true <- 1
beta1.true <- 1
yp <- c(1, 1, 0, 0)
xp <- c(1, 0, 1, 0)
# P(Y=y,X=x) in order of 11, 10, 01, 00
pi.true <- exp(yp * (beta0.true + beta1.true * xp))/(1 + exp(beta0.true + beta1.true * xp)) * p1.true^xp * (1 - p1.true)^(1 - xp)

## function that returns D values that produce max and min \delta
DD <- function(para) {
  eps <- 0.001
  SNY1 <- para[1]
  SNY0 <- para[2]
  SPY1 <- para[3]
  SPY0 <- para[4]
  SNX <- para[5]
  SPX <- para[6]
  L1 <- max(-SNY1 * SNX, -(1 - SNY1) * (1 - SNX))
  U1 <- min((1 - SNY1) * SNX, SNY1 * (1 - SNX))
  L2 <- max(-SNY0 * SPX, -(1 - SNY0) * (1 - SPX))
  U2 <- min((1 - SNY0) * SPX, SNY0 * (1 - SPX))
  L3 <- max(-SPY1 * SNX, -(1 - SPY1) * (1 - SNX))
  U3 <- min((1 - SPY1) * SNX, SPY1 * (1 - SNX))
  L4 <- max(-SPY0 * SPX, -(1 - SPY0) * (1 - SPX))
  U4 <- min((1 - SPY0) * SPX, SPY0 * (1 - SPX))
  obj <- function(D) {
    ans <- abs(sum(pi.true * D * c(1, -1, -1, 1)))
  }
  obj_minus <- function(D) {
    ans <- -abs(sum(pi.true * D * c(1, -1, -1, 1)))
  }
  Dlow <- optim(c(0, 0.01, 0.01, 1), obj, lower = c(L1, L2, L3, L4) + eps, upper = c(U1, U2, U3, U4) - eps, method = "L-BFGS-B")
  Dhigh <- optim(c(0, 0.01, 0.01, 1), obj_minus, lower = c(L1, L2, L3, L4) + eps, upper = c(U1, U2, U3, U4) - eps, method = "L-BFGS-B")
  ans <- list(D_l = round(c(Dlow$par, Dlow$val), 3), D_h = round(c(Dhigh$par, -Dhigh$val), 3))
  print(ans)
}

#### columns in the order of SNY1, SNY0, SPY1, SPY0, SNX, SPX
SNSP <- matrix(c(0.9, 0.85, 0.9, 0.85, 0.9, 0.9, 0.85, 0.8, 0.85, 0.8, 0.85, 0.85, 0.75, 0.7, 0.75, 0.7, 0.75, 0.75), nrow = 3, byrow = TRUE)

allcases <- matrix(NA, nrow = 6, ncol = 10)

for (s in 1:3) {
  DDtmp <- DD(SNSP[s, ])
  allcases[(s - 1) * 2 + 1, ] <- c(SNSP[s, ], DDtmp$D_l[1:4])
  allcases[2 * s, ] <- c(SNSP[s, ], DDtmp$D_h[1:4])
}

### data generation
Nv = 1000
Nm = 9000

para <- allcases[2, ]  # SN high,SP high,D low

SNY1 <- SPY1 <- SNX <- SPX <- para[1]
SNY0 <- SPY0 <- para[2]
SNY <- c(SNY1, SNY0)
SPY <- c(SPY1, SPY0)

DD <- para[-(1:6)]

C1 <- C0 <- Dsign <- matrix(NA, nrow = 4, ncol = 4)
for (i in 1:4) {
  # i index observed values j index true values
  for (j in 1:4) {
    ystar <- 2 - ceiling(i/2)
    xstar <- ifelse(floor(i/2) == ceiling(i/2), 0, 1)
    y <- 2 - ceiling(j/2)
    x <- ifelse(floor(j/2) == ceiling(j/2), 0, 1)
    
    C0[i, j] <- (SNY[2 - x]^ystar * (1 - SNY[2 - x])^(1 - ystar))^y * ((1 - SPY[2 - x])^ystar * SPY[2 - x]^(1 - ystar))^(1 - y) * (SNX^xstar * (1 -
                                                                                                                                                  SNX)^(1 - xstar))^x * ((1 - SPX)^xstar * SPX^(1 - xstar))^(1 - x)
    
    Dsign[i, j] <- (-1)^((x == xstar) + (y == ystar))
    
    C1[i, j] <- C0[i, j] + DD[j] * Dsign[i, j]
  }  # end of j loop
}  # end of 

C1 <- ifelse(C1 < 1e-15, 0, C1)
# split into validation datasets and main datasets

# validation data (ystar,xstar|y,x)(y,x) by columns (y,x) in order of 11, 10, 01, 00
pjoint.true <- as.vector(t(t(C1) * pi.true))
pstar.true <- C1 %*% pi.true[1:4]

set.seed(123)
validation <- as.vector(rmultinom(1, Nv, pjoint.true))
main <- as.vector(rmultinom(1, Nm, pstar.true))

xlabel <- ylabel <- 1:0
all <- as.data.frame(expand.grid(xlabel, ylabel, xlabel, ylabel))

names(all) <- c("xstar", "ystar", "x", "y")
all$count <- validation
SNY0_v <- sum(all$count[all$ystar == 1 & all$y == 1 & all$x == 0])/sum(all$count[all$y == 1 & all$x == 0])
SPY0_v <- sum(all$count[all$ystar == 0 & all$y == 0 & all$x == 0])/sum(all$count[all$y == 0 & all$x == 0])
SNY1_v <- sum(all$count[all$ystar == 1 & all$y == 1 & all$x == 1])/sum(all$count[all$y == 1 & all$x == 1])
SPY1_v <- sum(all$count[all$ystar == 0 & all$y == 0 & all$x == 1])/sum(all$count[all$y == 0 & all$x == 1])

SNX_v <- sum(all$count[all$xstar == 1 & all$x == 1])/sum(all$count[all$x == 1])
SPX_v <- sum(all$count[all$xstar == 0 & all$x == 0])/sum(all$count[all$x == 0])

D0 <- rep(NA, 1)
D0[1] <- all$count[all$ystar == 1 & all$xstar == 1 & all$y == 1 & all$x == 1]/sum(all$count[all$y == 1 & all$x == 1]) - SNY1_v * SNX_v
D0[2] <- all$count[all$ystar == 1 & all$xstar == 0 & all$y == 1 & all$x == 0]/sum(all$count[all$y == 1 & all$x == 0]) - SNY0_v * SPX_v
D0[3] <- all$count[all$ystar == 0 & all$xstar == 1 & all$y == 0 & all$x == 1]/sum(all$count[all$y == 0 & all$x == 1]) - SPY1_v * SNX_v
D0[4] <- all$count[all$ystar == 0 & all$xstar == 0 & all$y == 0 & all$x == 0]/sum(all$count[all$y == 0 & all$x == 0]) - SPY0_v * SPX_v
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
data <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, a1 = a1, b1 = b1)

########### differential Y, nondifferential X, and dependent model

ddmodelx <- "
model{
SNY <- c(SNY1, SNY0)
SPY <- c(SPY1, SPY0)
for (i in 1:4) {
    # i index observed values j index true values
    for (j in 1:4) {
        ystar[i, j] <- 2 - round(i/2)
        xstar[i, j] <- ifelse(round(i/2) == trunc(i/2), 0, 1)
        y[i, j] <- 2 - round(j/2)
        x[i, j] <- ifelse(round(j/2) == trunc(j/2), 0, 1)
        # C[i,j]=P(Ystar=ystar,Xstar=xstar|Y=y,X=x); i indexing the pair (ystar,xstar) and j indexing the pair (y,x)
        C[i, j] <- (SNY[2 - x[i, j]]^ystar[i, j] * (1 - SNY[2 - x[i, j]])^(1 - ystar[i, j]))^y[i, j] * ((1 - SPY[2 - x[i, j]])^ystar[i, j] * SPY[2 - x[i,
            j]]^(1 - ystar[i, j]))^(1 - y[i, j]) * (SNX^xstar[i, j] * (1 - SNX)^(1 - xstar[i, j]))^x[i, j] * ((1 - SPX)^xstar[i, j] * SPX^(1 - xstar[i,
            j]))^(1 - x[i, j]) + D[j] * (-1)^(ifelse(ystar[i, j] == y[i, j], 1, 0) + ifelse(xstar[i, j] == x[i, j], 1, 0))
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
SNY1 ~ dunif(1 - SPY1, 1)
SPY1 ~ dunif(0.5, 1)

SNY0 ~ dunif(1 - SPY0, 1)
SPY0 ~ dunif(0.5, 1)

L[1] <- max(-SNY1 * SNX, -(1 - SNY1) * (1 - SNX))
U[1] <- min((1 - SNY1) * SNX, SNY1 * (1 - SNX))
L[2] <- max(-SNY0 * SPX, -(1 - SNY0) * (1 - SPX))
U[2] <- min((1 - SNY0) * SPX, SNY0 * (1 - SPX))
L[3] <- max(-SPY1 * SNX, -(1 - SPY1) * (1 - SNX))
U[3] <- min((1 - SPY1) * SNX, SPY1 * (1 - SNX))
L[4] <- max(-SPY0 * SPX, -(1 - SPY0) * (1 - SPX))
U[4] <- min((1 - SPY0) * SPX, SPY0 * (1 - SPX))

# priors for D parameters
for (j in 1:4) {
    Dt[j] ~ dunif(0, 1)
    D[j] <- L[j] + (U[j] - L[j]) * Dt[j]
 }    
}#end of model
"

monitor <- c("SNY", "SNX", "SPY", "SPX", "D", "beta0", "beta1")  # names of all the parameters
inits1 <- list(SNY0 = SNY0_v, SPY0 = SPY0_v, SNY1 = SNY1_v, SPY1 = SPY1_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, Dt = rep(0.8, 4), beta0 = beta0_ini, beta1 = beta1_ini,
               .RNG.name = "base::Super-Duper", .RNG.seed = 1)
inits2 <- list(SNY0 = SNY0_v, SPY0 = SPY0_v, SNY1 = SNY1_v, SPY1 = SPY1_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, Dt = rep(0.5, 4), beta0 = beta0_ini, beta1 = beta1_ini,
               .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2)

opt_dd <- run.jags(model = ddmodelx, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000, method = "rjags", inits = list(inits1,
                                                                                                                                                               inits2))

out_dd <- summary(opt_dd)[, c(2, 4, 5, 11)]

###########independent model

dimodel <- "
model{
SNY <- c(SNY1, SNY0)
SPY <- c(SPY1, SPY0)

for (i in 1:4) {
    # i index observed values j index true values
    for (j in 1:4) {
        ystar[i, j] <- 2 - round(i/2)
        xstar[i, j] <- ifelse(round(i/2) == trunc(i/2), 0, 1)
        y[i, j] <- 2 - round(j/2)
        x[i, j] <- ifelse(round(j/2) == trunc(j/2), 0, 1)
        # C[i,j]=P(Ystar=ystar,Xstar=xstar|Y=y,X=x); i indexing the pair (ystar,xstar) and j indexing the pair (y,x)
        C[i, j] <- (SNY[2 - x[i, j]]^ystar[i, j] * (1 - SNY[2 - x[i, j]])^(1 - ystar[i, j]))^y[i, j] * ((1 - SPY[2 - x[i, j]])^ystar[i, j] * SPY[2 - x[i,
            j]]^(1 - ystar[i, j]))^(1 - y[i, j]) * (SNX^xstar[i, j] * (1 - SNX)^(1 - xstar[i, j]))^x[i, j] * ((1 - SPX)^xstar[i, j] * SPX^(1 - xstar[i,
            j]))^(1 - x[i, j])
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

SNY1 ~ dunif(1 - SPY1, 1)
SPY1 ~ dunif(0.5, 1)
SNY0 ~ dunif(1 - SPY0, 1)
SPY0 ~ dunif(0.5, 1)    
}
"

monitor <- c("SNY", "SNX", "SPY", "SPX", "beta0", "beta1")  # names of all the parameters
# note both p and cov are in the ordering of 11,10,01,00

inits1 <- list(SNY0 = SNY0_v, SPY0 = SPY0_v, SNY1 = SNY1_v, SPY1 = SPY1_v, SNX = SNX_v, SPX = SPX_v, p1 = px0, beta0 = beta0_ini, beta1 = beta1_ini, .RNG.name = "base::Super-Duper",
               .RNG.seed = 1)
inits2 <- list(SNY0 = SNY0_v, SPY0 = SPY0_v, SNY1 = SNY1_v, SPY1 = SPY1_v, SNX = SNX_v, SPX = SPX_v, p1 = 0.2, beta0 = beta0_ini, beta1 = beta1_ini, .RNG.name = "base::Wichmann-Hill",
               .RNG.seed = 2)

opt_di <- run.jags(model = dimodel, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000, method = "rjags", inits = list(inits1,inits2))
                                                                                                                                                              

out_di <- summary(opt_di)[, c(2, 4, 5, 11)]  # what does mode mean

############################### naive model

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
        C[i, j] <- exp(ystar[i, j] * (beta0 + xstar[i, j] * beta1))/(1 + exp(beta0 + xstar[i, j] * beta1)) * p1^(xstar[i, j]) * (1 - p1)^(1 - xstar[i,
            j])
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

monitor <- c("beta0", "beta1")  # names of all the parameters
# note both p and cov are in the ordering of 11,10,01,00

inits1 <- list(beta0 = beta0_ini, beta1 = beta1_ini, p1 = px0, .RNG.name = "base::Super-Duper", .RNG.seed = 1)
inits2 <- list(beta0 = 0, beta1 = 1, p1 = 0.2, .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2)

opt_naive <- run.jags(model = naive, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000, burnin = 5000, method = "rjags", inits = list(inits1,
                                                                                                                                                               inits2))

out_naive <- summary(opt_naive)[, c(2, 4, 5, 11)]