# categorical/ordinal Y and X, each with three categories
# both subject to nondifferential misclassification errors
# misclassification only occurs between two adjacent categories

### data generation
######################################
Nv <- 1000
Nm <- 9000
CY <- 0.9
CX <- 0.9
alpha <- c(-1, 1)
beta <- c(0, 1, 2)
px <- c(1 / 3, 1 / 3, 1 / 3) # pmf of X

I <- 3
J <- 3
xlabel <- 0:(I - 1)
ylabel <- 0:(J - 1)
all <- as.data.frame(expand.grid(ylabel, xlabel, ylabel, xlabel))
names(all) <- c("x", "y", "xstar", "ystar")
all_order <- all[, 4:1]

yx <- as.data.frame(expand.grid(ylabel, xlabel))
names(yx) <- c("x", "y")
yx <- data.matrix(yx)

expit <- function(t) {
  exp(t) / (1 + exp(t))
}

prob_vec <- function(x, y) { # x and y are vectors labeld 0:2
  cond <- joint <- rep(0, length(x)) ## cond P(Y|X)
  for (i in 1:length(x)) {
    if (y[i] == 0) {
      cond[i] <- expit(alpha[1] + beta[x[i] + 1])
    } else if (y[i] == I - 1) {
      cond[i] <- 1 - expit(alpha[I - 1] + beta[x[i] + 1])
    } else {
      cond[i] <- expit(alpha[y[i] + 1] + beta[x[i] + 1]) - expit(alpha[y[i]] + beta[x[i] + 1])
    } # end of else
  } # end of i loop
  joint <- cond * px[x + 1]
  return(joint)
}

# generate Y and X first
pi.true <- prob_vec(yx[, 1], yx[, 2])

mis_matrix <- function(pcorrect, dim) {
  out <- diag(pcorrect, nrow = dim)
  out[1:2, 1] <- c(pcorrect, 1 - pcorrect)
  out[(dim - 1):dim, dim] <- c(1 - pcorrect, pcorrect)
  for (j in 2:(dim - 1)) {
    out[j - 1, j] <- out[j + 1, j] <- (1 - pcorrect) / 2
  }
  return(out)
}

MX <- mis_matrix(CX, I)
MY <- mis_matrix(CY, J)
total <- I * J

C1 <- C0 <- matrix(NA, nrow = total, ncol = total)
# p*=C%*%p
# C1 is the multiplication matrix under H_1
# C0 is the multiplication matrix under H_0


# 16 free D parameters in D matrix
# output matrix D^2/(MX*MY), elementwise operation

Dmatrix <- function(dvec, CY, CX) { # function to convert the 16 free D parameters to matrix
  MY <- mis_matrix(CY, I)
  MX <- mis_matrix(CX, J)
  D <- phi <- matrix(NA, nrow = total, ncol = total)
  for (k in 1:total) {
    for (m in 1:total) {
      ystar <- floor((k - 1) / I)
      xstar <- k - 1 - ystar * J
      y <- floor((m - 1) / J)
      x <- m - 1 - y * J
      if (abs(ystar - y) > 1 | abs(xstar - x) > 1) {
        D[k, m] <- phi[k, m] <- 0
      } else {
        phi[k, m] <- 1 / (MY[ystar + 1, y + 1] * MX[xstar + 1, x + 1])
      }
    }
  }
  D[is.na(D[, 1]), 1] <- c(dvec[1], -dvec[1], -dvec[1], dvec[1])
  D[is.na(D[, 2]), 2] <- c(dvec[2:3], -sum(dvec[2:3]), -dvec[2:3], sum(dvec[2:3]))
  D[is.na(D[, 3]), 3] <- c(dvec[4], -dvec[4], -dvec[4], dvec[4])
  D[is.na(D[, 4]), 4] <- c(dvec[5], -dvec[5], dvec[6], -dvec[6], -sum(dvec[5:6]), sum(dvec[5:6]))
  D[is.na(D[, 5]), 5] <- c(dvec[7:8], -dvec[7] - dvec[8], dvec[9], dvec[10], -sum(dvec[9:10]), -sum(dvec[c(7, 9)]), -sum(dvec[c(8, 10)]), sum(dvec[7:10]))
  D[is.na(D[, 6]), 6] <- c(dvec[11], -dvec[11], dvec[12], -dvec[12], -sum(dvec[11:12]), sum(dvec[11:12]))
  D[is.na(D[, 7]), 7] <- c(dvec[13], -dvec[13], -dvec[13], dvec[13])
  D[is.na(D[, 8]), 8] <- c(dvec[14:15], -sum(dvec[14:15]), -dvec[14:15], sum(dvec[14:15]))
  D[is.na(D[, 9]), 9] <- c(dvec[16], -dvec[16], -dvec[16], dvec[16])

  Dphi <- D^2 * phi # squared difference
  ans <- (list(D = D, Dphi = Dphi))
}

## low and high dependence defined by normalized Cramer's V, i.e., phi_V

## optimization subject to constraints

library(alabama)
## function to get the 16 D parameters so as to have high dependence
DD <- function(para) {
  CY <- para[1]
  CX <- para[2]
  MX <- mis_matrix(CX, J)
  MY <- mis_matrix(CY, I)

  # upper and lower limites of 16 D parameters
  boundary <- function(y, x, ystar, xstar) {
    L <- max(-(1 - MY[ystar + 1, y + 1]) * (1 - MX[xstar + 1, x + 1]), -MY[ystar + 1, y + 1] * MX[xstar + 1, x + 1])
    U <- min(MY[ystar + 1, y + 1] * (1 - MX[xstar + 1, x + 1]), (1 - MY[ystar + 1, y + 1]) * MX[xstar + 1, x + 1])
    return(c(L, U))
  }
  ## boundaries for d1
  r1 <- boundary(0, 0, 0, 0)
  r2 <- -1 * (boundary(0, 0, 0, 1))[2:1]
  r3 <- -1 * (boundary(0, 0, 1, 0))[2:1]
  r4 <- (boundary(0, 0, 1, 1))
  rm <- cbind(r1, r2, r3, r4)
  L1 <- max(rm[1, ])
  U1 <- min(rm[2, ])

  # boundaries for d2
  r1 <- boundary(0, 1, 0, 0)
  r2 <- -1 * (boundary(0, 1, 1, 0))[2:1]
  rm <- cbind(r1, r2)
  L2 <- max(rm[1, ])
  U2 <- min(rm[2, ])

  # boundary for d3
  r1 <- boundary(0, 1, 0, 1)
  r2 <- -1 * (boundary(0, 1, 1, 1))[2:1]
  rm <- cbind(r1, r2)
  L3 <- max(rm[1, ])
  U3 <- min(rm[2, ])
  # boundary for d2+d3
  s1 <- -1 * (boundary(0, 1, 0, 2))[2:1]
  s2 <- boundary(0, 1, 1, 2)
  sm <- cbind(s1, s2)
  LS12 <- max(sm[1, ])
  US12 <- min(sm[2, ])

  # boundary for d4
  r1 <- boundary(0, 2, 0, 1)
  r2 <- -1 * (boundary(0, 2, 0, 2))[2:1]
  r3 <- -1 * (boundary(0, 2, 1, 1))[2:1]
  r4 <- boundary(0, 2, 1, 2)
  rm <- cbind(r1, r2, r3, r4)
  L4 <- max(rm[1, ])
  U4 <- min(rm[2, ])

  # boundary for d5
  r1 <- boundary(1, 0, 0, 0)
  r2 <- -1 * (boundary(1, 0, 0, 1))[2:1]
  rm <- cbind(r1, r2)
  L5 <- max(rm[1, ])
  U5 <- min(rm[2, ])

  # boundary for d6
  r1 <- boundary(1, 0, 1, 0)
  r2 <- -1 * (boundary(1, 0, 1, 1))[2:1]
  rm <- cbind(r1, r2)
  L6 <- max(rm[1, ])
  U6 <- min(rm[2, ])

  # boundary for d5+d6
  s1 <- -1 * (boundary(1, 0, 2, 0))[2:1]
  s2 <- boundary(1, 0, 2, 1)
  sm <- cbind(s1, s2)
  LS56 <- max(sm[1, ])
  US56 <- min(sm[2, ])

  # boundary for d7
  r1 <- boundary(1, 1, 0, 0)
  L7 <- r1[1]
  U7 <- r1[2]

  # boundary for d8
  r1 <- boundary(1, 1, 0, 1)
  L8 <- r1[1]
  U8 <- r1[2]

  # boundary for d9
  r1 <- boundary(1, 1, 1, 0)
  L9 <- r1[1]
  U9 <- r1[2]

  # boundary for d10
  r1 <- boundary(1, 1, 1, 1)
  L10 <- r1[1]
  U10 <- r1[2]

  ######## update boundaries for d7-10
  r78 <- -1 * boundary(1, 1, 0, 2)[2:1]
  r910 <- -1 * boundary(1, 1, 1, 2)[2:1]
  r79 <- -1 * boundary(1, 1, 2, 0)[2:1]
  r810 <- -1 * boundary(1, 1, 2, 1)[2:1]
  s79810 <- boundary(1, 1, 2, 2)

  # boundary for d11
  r1 <- boundary(1, 2, 0, 1)
  r2 <- -1 * (boundary(1, 2, 0, 2))[2:1]
  rm <- cbind(r1, r2)
  L11 <- max(rm[1, ])
  U11 <- min(rm[2, ])

  # boundary for d12
  r1 <- boundary(1, 2, 1, 1)
  r2 <- -1 * (boundary(1, 2, 1, 2))[2:1]
  rm <- cbind(r1, r2)
  L12 <- max(rm[1, ])
  U12 <- min(rm[2, ])
  # boundary for d11+d12
  s1 <- -1 * (boundary(1, 2, 2, 1))[2:1]
  s2 <- boundary(1, 2, 2, 2)
  sm <- cbind(s1, s2)
  LS1112 <- max(sm[1, ])
  US1112 <- min(sm[2, ])

  # boundary for d13
  r1 <- boundary(2, 0, 1, 0)
  r2 <- -1 * (boundary(2, 0, 1, 1))[2:1]
  r4 <- boundary(2, 0, 2, 1)
  r3 <- -1 * (boundary(2, 0, 2, 0))[2:1]
  rm <- cbind(r1, r2, r3, r4)
  L13 <- max(rm[1, ])
  U13 <- min(rm[2, ])

  # boundary for d14
  r1 <- boundary(2, 1, 1, 0)
  r2 <- -1 * (boundary(2, 1, 2, 0))[2:1]
  rm <- cbind(r1, r2)
  L14 <- max(rm[1, ])
  U14 <- min(rm[2, ])

  # boundary for d15
  r1 <- boundary(2, 1, 1, 1)
  r2 <- -1 * (boundary(2, 1, 2, 1))[2:1]
  rm <- cbind(r1, r2)
  L15 <- max(rm[1, ])
  U15 <- min(rm[2, ])
  # boundary for d14+d15
  s1 <- -1 * (boundary(2, 1, 1, 2))[2:1]
  s2 <- boundary(2, 1, 2, 2)
  sm <- cbind(s1, s2)
  LS1415 <- max(sm[1, ])
  US1415 <- min(sm[2, ])

  # boundary for d16
  r1 <- boundary(2, 2, 1, 1)
  r2 <- -1 * (boundary(2, 2, 1, 2))[2:1]
  r4 <- boundary(2, 2, 2, 2)
  r3 <- -1 * (boundary(2, 2, 2, 1))[2:1]
  rm <- cbind(r1, r2, r3, r4)
  L16 <- max(rm[1, ])
  U16 <- min(rm[2, ])

  DDM <- cbind(
    c(L1, U1), c(L2, U2), c(L3, U3), c(L4, U4), c(L5, U5), c(L6, U6), c(L7, U7), c(L8, U8), c(L9, U9), c(L10, U10),
    c(L11, U11), c(L12, U12), c(L13, U13), c(L14, U14), c(L15, U15), c(L16, U16)
  )
  ## objective function associated with maximum phi_V. Note maximization is minimizing minus phi_V
  objh <- function(d) {
    Phi <- Dmatrix(d, CY, CX)$Dphi
    ans <- -1 * sum(pi.true * sqrt(apply(Phi, 2, sum)) / sqrt(2))
  }

  ## objective function associated with minimum phi_V. Note maximization is minimizing minus phi_V
  objl <- function(d) {
    Phi <- Dmatrix(d, CY, CX)$Dphi
    ans <- sum(pi.true * sqrt(apply(Phi, 2, sum)) / sqrt(2))
  }

  eps <- 0.001
  hin <- function(d) {
    h <- rep(NA, 1)
    for (k in 1:16) {
      h[(k - 1) * 2 + 1] <- d[k] - (DDM[1, k] + eps)
      h[(k - 1) * 2 + 2] <- DDM[2, k] - eps - d[k]
    }
    h[33] <- d[2] + d[3] - (LS12 + eps)
    h[34] <- US12 - eps - d[2] - d[3]
    h[35] <- d[5] + d[6] - (LS56 + eps)
    h[36] <- US56 - eps - d[5] - d[6]
    h[37] <- d[11] + d[12] - (LS1112 + eps)
    h[38] <- US1112 - eps - d[11] - d[12]
    h[39] <- d[14] + d[15] - (LS1415 + eps)
    h[40] <- US1415 - eps - d[14] - d[15]
    h[41] <- d[7] + d[8] - (r78[1] + eps)
    h[42] <- r78[2] - eps - d[7] - d[8]
    h[43] <- d[9] + d[10] - (r910[1] + eps)
    h[44] <- r910[2] - eps - d[9] - d[10]
    h[45] <- d[7] + d[9] - (r79[1] + eps)
    h[46] <- r79[2] - eps - d[7] - d[9]
    h[47] <- d[8] + d[10] - (r810[1] + eps)
    h[48] <- r810[2] - eps - d[8] - d[10]
    h[49] <- d[7] + d[8] + d[9] + d[10] - (s79810[1] + eps)
    h[50] <- s79810[2] - eps - d[7] - d[8] - d[9] - d[10]
    h
  }

  Dhigh <- auglag(par = rep(0.001, 16), fn = objh, hin = hin, control.outer = list(trace = FALSE))
  Dlow <- auglag(par = rep(0.001, 16), fn = objl, hin = hin, control.outer = list(trace = FALSE))
  ans <- list(D_h = round(c(Dhigh$par, -Dhigh$value), 4), D_l = round(c(Dlow$par, Dlow$value), 4))
  print(ans)
}

### generate validation data
out <- DD(c(CY, CX))
dh <- out$D_h[1:16]
phiDh <- out$D_h[17]

dl <- out$D_l[1:16]
phiDl <- out$D_l[17]


Dout <- Dmatrix(dh, CY, CX)
D <- Dout$D
for (k in 1:total) {
  for (m in 1:total) {
    ystar <- floor((k - 1) / I)
    xstar <- k - 1 - ystar * J
    y <- floor((m - 1) / J)
    x <- m - 1 - y * J
    print(c(ystar, xstar, y, x))
    C0[k, m] <- MY[ystar + 1, y + 1] * MX[xstar + 1, x + 1]
    C1[k, m] <- C0[k, m] + D[k, m]
  } # end of m loop
} # end of k loop


A1 <- MY %x% MX + D

## (y,x) in order of 00,01,02,10,11,12,20,21,22
## (ystar,xstar|y,x) in order of 00,01,02,10,11,12,20,21,22 (ystar, xstar) for row and (y,x) for column
pstar_true <- A1 %*% pi.true
pjoint_true <- as.vector(t(A1) * pi.true)

set.seed(123)
validation <- as.vector(rmultinom(1, Nv, pjoint_true))
main <- as.vector(rmultinom(1, Nm, pstar_true))

library(coda)
library(runjags)

### data
px0 <- px
a1 <- 116
b1 <- 12
data <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, px0 = px0, a1 = a1, b1 = b1)

ddmodel <- "
model{
## true joint for (Y,X) in order of 00,01,02,10,11,12,20,21,22
p[1] <- 1 / (1 + exp(-alpha1)) * px[1]
p[2] <- 1 / (1 + exp(-alpha1 - beta1)) * px[2]
p[3] <- 1 / (1 + exp(-alpha1 - beta2)) * px[3]
p[4] <- (1 / (1 + exp(-alpha2)) - 1 / (1 + exp(-alpha1))) * px[1]
p[5] <- (1 / (1 + exp(-alpha2 - beta1)) - 1 / (1 + exp(-alpha1 - beta1))) * px[2]
p[6] <- (1 / (1 + exp(-alpha2 - beta2)) - 1 / (1 + exp(-alpha1 - beta2))) * px[3]
p[7] <- (1 - 1 / (1 + exp(-alpha2))) * px[1]
p[8] <- (1 - 1 / (1 + exp(-alpha2 - beta1))) * px[2]
p[9] <- (1 - 1 / (1 + exp(-alpha2 - beta2))) * px[3]

# MY
MY[1:3, 1] <- c(CY, 1 - CY, 0)
MY[1:3, 2] <- c((1 - CY) / 2, CY, (1 - CY) / 2)
MY[1:3, 3] <- c(0, 1 - CY, CY)
# MX
MX[1:3, 1] <- c(CX, 1 - CX, 0)
MX[1:3, 2] <- c((1 - CX) / 2, CX, (1 - CX) / 2)
MX[1:3, 3] <- c(0, 1 - CX, CX)

# D matrix
D[1:9, 1] <- c(dvec[1], -dvec[1], 0, -dvec[1], dvec[1], 0, 0, 0, 0)
D[1:9, 2] <- c(dvec[2:3], -sum(dvec[2:3]), -dvec[2:3], sum(dvec[2:3]), 0, 0, 0)
D[1:9, 3] <- c(0, dvec[4], -dvec[4], 0, -dvec[4], dvec[4], 0, 0, 0)
D[1:9, 4] <- c(dvec[5], -dvec[5], 0, dvec[6], -dvec[6], 0, -sum(dvec[5:6]), sum(dvec[5:6]), 0)
D[1:9, 5] <- c(dvec[7:8], -dvec[7] - dvec[8], dvec[9], dvec[10], -sum(dvec[9:10]), -sum(dvec[c(7, 9)]), -sum(dvec[c(8, 10)]), sum(dvec[7:10]))
D[1:9, 6] <- c(0, dvec[11], -dvec[11], 0, dvec[12], -dvec[12], 0, -sum(dvec[11:12]), sum(dvec[11:12]))
D[1:9, 7] <- c(0, 0, 0, dvec[13], -dvec[13], 0, -dvec[13], dvec[13], 0)
D[1:9, 8] <- c(0, 0, 0, dvec[14:15], -sum(dvec[14:15]), -dvec[14:15], sum(dvec[14:15]))
D[1:9, 9] <- c(0, 0, 0, 0, dvec[16], -dvec[16], 0, -dvec[16], dvec[16])

for (k in 1:9) {
  for (m in 1:9) {
    ystar[k, m] <- trunc((k - 1) / 3)
    xstar[k, m] <- k - 1 - ystar[k, m] * 3
    y[k, m] <- trunc((m - 1) / 3)
    x[k, m] <- m - 1 - y[k, m] * 3
    QY[k, m] <- MY[ystar[k, m] + 1, y[k, m] + 1]
    QX[k, m] <- MX[xstar[k, m] + 1, x[k, m] + 1]
    C1[k, m] <- QY[k, m] * QX[k, m] + D[k, m]
    L[k, m] <- max(-(1 - QY[k, m]) * (1 - QX[k, m]), -QY[k, m] * QX[k, m])
    U[k, m] <- min(QY[k, m] * (1 - QX[k, m]), QX[k, m] * (1 - QY[k, m]))
  }
}

# boundaries for the 16 D
DL[1] <- max(L[1, 1], -U[2, 1], -U[4, 1], L[5, 1])
DU[1] <- min(U[1, 1], -L[2, 1], -L[4, 1], U[5, 1])

DL[2] <- max(L[1, 2], -U[4, 2])
DU[2] <- min(U[1, 2], -L[4, 2])

DL[3] <- max(L[2, 2], -U[5, 2])
DU[3] <- min(U[2, 2], -L[5, 2])

DL[4] <- max(L[2, 3], -U[3, 3], -U[5, 3], L[6, 3])
DU[4] <- min(U[2, 3], -L[3, 3], -L[5, 3], U[6, 3])

DL[5] <- max(L[1, 4], -U[2, 4])
DU[5] <- min(U[1, 4], -L[2, 4])

DL[6] <- max(L[4, 4], -U[5, 4])
DU[6] <- min(U[4, 4], -L[5, 4])

DL[7] <- L[1, 5]
DU[7] <- U[1, 5]

DL[8] <- L[2, 5]
DU[8] <- U[2, 5]

DL[9] <- L[4, 5]
DU[9] <- U[4, 5]

DL[10] <- L[5, 5]
DU[10] <- U[5, 5]

DL[11] <- max(L[2, 6], -U[3, 6])
DU[11] <- min(U[2, 6], -L[3, 6])

DL[12] <- max(L[5, 6], -U[6, 6])
DU[12] <- min(U[5, 6], -L[6, 6])

DL[13] <- max(L[4, 7], -U[5, 7])
DU[13] <- min(U[4, 7], -L[5, 7])

DL[14] <- max(L[4, 8], -U[7, 8])
DU[14] <- min(U[4, 8], -L[7, 8])

DL[15] <- max(L[5, 8], -U[8, 8])
DU[15] <- min(U[5, 8], -L[8, 8])

DL[16] <- max(L[5, 9], -U[6, 9], -U[8, 9], L[9, 9])
DU[16] <- min(U[5, 9], -L[6, 9], -L[8, 9], U[9, 9])

# validation data
validation[1:81] ~ dmulti(pjoint[1:81], Nv)

for (s in 1:9) {
  pjoint[(9 * (s - 1) + 1):(9 * s)] <- C1[s, ] * p
}

# main data
main[1:9] ~ dmulti(pstar[1:9], Nm)
pstar[1:9] <- C1[1:9, 1:9] %*% p[1:9]

# priors
CX ~ dbeta(a1, b1)
CY ~ dbeta(a1, b1)

beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)
alpha1 ~ dnorm(0, 0.001)T(, alpha2)
alpha2 ~ dnorm(0, 0.01)

px[1:3] ~ ddirch(kappa * px0[1:3]) # px0 is estimated marginal distribution of x based on the validation data
kappa ~ dgamma(2, 0.001)

for (j in 1:16) {
  dvec[j] ~ dunif(DL[j], DU[j])
}

# constraints
ind23 <- step(sum(dvec[2:3]) - L[6, 2]) * step(U[6, 2] - sum(dvec[2:3])) * step(-sum(dvec[2:3]) - L[3, 2]) * step(U[3, 2] + sum(dvec[2:3]))
ind56 <- step(sum(dvec[5:6]) - L[8, 4]) * step(U[8, 4] - sum(dvec[5:6])) * step(-sum(dvec[5:6]) - L[7, 4]) * step(U[7, 4] + sum(dvec[5:6]))
ind78 <- step(-sum(dvec[7:8]) - L[3, 5]) * step(U[3, 5] + sum(dvec[7:8]))
ind910 <- step(-sum(dvec[9:10]) - L[6, 5]) * step(U[6, 5] + sum(dvec[9:10]))
ind79 <- step(-dvec[7] - dvec[9] - L[7, 5]) * step(U[7, 5] + dvec[7] + dvec[9])
ind810 <- step(-dvec[8] - dvec[10] - L[8, 5]) * step(U[8, 5] + dvec[8] + dvec[10])
ind78910 <- step(sum(dvec[7:10]) - L[9, 5]) * step(U[9, 5] - sum(dvec[7:10]))
ind1112 <- step(sum(dvec[11:12]) - L[9, 6]) * step(U[9, 6] - sum(dvec[11:12])) * step(-sum(dvec[11:12]) - L[8, 6]) * step(U[8, 6] + sum(dvec[11:12]))
ind1415 <- step(sum(dvec[14:15]) - L[9, 8]) * step(U[9, 8] - sum(dvec[14:15])) * step(-sum(dvec[14:15]) - L[6, 8]) * step(U[6, 8] + sum(dvec[14:15]))
ind <- ind23 * ind56 * ind78 * ind910 * ind79 * ind810 * ind78910 * ind1112 * ind1415
}
"

monitor <- c("CY", "CX", "dvec", "beta1", "beta2", "alpha1", "alpha2", "ind") # names of all the parameters
inits1 <- list(
  CY = 0.8, CX = 0.8, beta1 = 1, beta2 = 1, alpha1 = -0.5, alpha2 = 1, dvec = rep(0, 16), kappa = 1,
  .RNG.name = "base::Super-Duper", .RNG.seed = 1
  )
inits2 <- list(
  CY = 0.8, CX = 0.8, beta1 = 0, beta2 = 1, alpha1 = -1, alpha2 = 1, dvec = rep(0.001, 16), kappa = 1,
  .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2
)

opt_dd <- run.jags(
  model = ddmodel, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000,
  burnin = 5000, method = "rjags", inits = list(inits1, inits2)
)

## save output
opt <- as.matrix(opt_dd$mcmc)
postsummary <- function(x) {
  ans <- c(median(x), mean(x), sd(x))
} ## median, mean and sd
out_dd <- t(apply(opt[opt[, 23] == 1, ], 2, postsummary)) ## only including samples that satisfy contraints for D parameters

##############################
# independent model
##############################
data <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, px0 = px0, a1 = a1, b1 = b1)

indmodel <- "
model{
## true joint for (Y,X) in order of 00,01,02,10,11,12,20,21,22
p[1] <- 1/(1 + exp(-alpha1)) * px[1]
p[2] <- 1/(1 + exp(-alpha1 - beta1)) * px[2]
p[3] <- 1/(1 + exp(-alpha1 - beta2)) * px[3]
p[4] <- (1/(1 + exp(-alpha2)) - 1/(1 + exp(-alpha1))) * px[1]
p[5] <- (1/(1 + exp(-alpha2 - beta1)) - 1/(1 + exp(-alpha1 - beta1))) * px[2]
p[6] <- (1/(1 + exp(-alpha2 - beta2)) - 1/(1 + exp(-alpha1 - beta2))) * px[3]
p[7] <- (1 - 1/(1 + exp(-alpha2))) * px[1]
p[8] <- (1 - 1/(1 + exp(-alpha2 - beta1))) * px[2]
p[9] <- (1 - 1/(1 + exp(-alpha2 - beta2))) * px[3]

# MY
MY[1:3, 1] <- c(CY, 1 - CY, 0)
MY[1:3, 2] <- c((1 - CY)/2, CY, (1 - CY)/2)
MY[1:3, 3] <- c(0, 1 - CY, CY)
# MX
MX[1:3, 1] <- c(CX, 1 - CX, 0)
MX[1:3, 2] <- c((1 - CX)/2, CX, (1 - CX)/2)
MX[1:3, 3] <- c(0, 1 - CX, CX)

for (k in 1:9) {
    for (m in 1:9) {
        ystar[k, m] = trunc((k - 1)/3)
        xstar[k, m] = k - 1 - ystar[k, m] * 3
        y[k, m] = trunc((m - 1)/3)
        x[k, m] = m - 1 - y[k, m] * 3
        QY[k, m] = MY[ystar[k, m] + 1, y[k, m] + 1]
        QX[k, m] = MX[xstar[k, m] + 1, x[k, m] + 1]
        C0[k, m] <- QY[k, m] * QX[k, m]
    }
}

# validation data
validation[1:81] ~ dmulti(pjoint[1:81], Nv)

for (s in 1:9) {
    pjoint[(9 * (s - 1) + 1):(9 * s)] <- C0[s, ] * p
}

# main data
main[1:9] ~ dmulti(pstar[1:9], Nm)
pstar[1:9] <- C0[1:9, 1:9] %*% p[1:9]

# priors
CX ~ dbeta(a1, b1)
CY ~ dbeta(a1, b1)

beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)
alpha1 ~ dnorm(0, 0.001)T(, alpha2)
alpha2 ~ dnorm(0, 0.01)

px[1:3] ~ ddirch(kappa * px0[1:3])  # px0 is estimated marginal distribution of x based on the validation data
kappa ~ dgamma(2, 0.001)
}
"
monitor <- c("CY", "CX", "beta1", "beta2", "alpha1", "alpha2") # names of all the parameters
inits1 <- list(
  CY = 0.8, CX = 0.8, beta1 = 1, beta2 = 1, alpha1 = -0.5, alpha2 = 1, kappa = 1,
  .RNG.name = "base::Super-Duper", .RNG.seed = 1
)
inits2 <- list(
  CY = 0.8, CX = 0.8, beta1 = 0, beta2 = 1, alpha1 = -1, alpha2 = 1, kappa = 1,
  .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2
)

opt_ind <- run.jags(
  model = indmodel, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000,
  burnin = 5000, method = "rjags", inits = list(inits1, inits2)
)

out_ind <- summary(opt_ind)[, c(2, 4, 5, 11)]


#########################
# naive model
##########################

data <- list(validation = validation, main = main, Nv = Nv, Nm = Nm, px0 = px0)

naivemodel <- "
model{
## true joint for (Y,X) in order of 00,01,02,10,11,12,20,21,22
p[1] <- 1/(1 + exp(-alpha1)) * px[1]
p[2] <- 1/(1 + exp(-alpha1 - beta1)) * px[2]
p[3] <- 1/(1 + exp(-alpha1 - beta2)) * px[3]
p[4] <- (1/(1 + exp(-alpha2)) - 1/(1 + exp(-alpha1))) * px[1]
p[5] <- (1/(1 + exp(-alpha2 - beta1)) - 1/(1 + exp(-alpha1 - beta1))) * px[2]
p[6] <- (1/(1 + exp(-alpha2 - beta2)) - 1/(1 + exp(-alpha1 - beta2))) * px[3]
p[7] <- (1 - 1/(1 + exp(-alpha2))) * px[1]
p[8] <- (1 - 1/(1 + exp(-alpha2 - beta1))) * px[2]
p[9] <- (1 - 1/(1 + exp(-alpha2 - beta2))) * px[3]

# validation data
validation[1:81] ~ dmulti(pjoint[1:81], Nv)

for (s in 1:9) {
    pjoint[(9 * (s - 1) + 1):(9 * s)] <- p[s] * p
}

# main data
main[1:9] ~ dmulti(p[1:9], Nm)

# priors
beta1 ~ dnorm(0, 0.001)
beta2 ~ dnorm(0, 0.001)
alpha1 ~ dnorm(0, 0.001)T(, alpha2)
alpha2 ~ dnorm(0, 0.01)

px[1:3] ~ ddirch(kappa * px0[1:3])  # px0 is estimated marginal distribution of x based on the validation data
kappa ~ dgamma(2, 0.001)
}
"
monitor <- c("beta1", "beta2", "alpha1", "alpha2") # names of all the parameters
inits1 <- list(
  beta1 = 1, beta2 = 1, alpha1 = -0.5, alpha2 = 1, kappa = 1,
  .RNG.name = "base::Super-Duper", .RNG.seed = 1
)
inits2 <- list(
  beta1 = 0, beta2 = 1, alpha1 = -1, alpha2 = 1, kappa = 1,
  .RNG.name = "base::Wichmann-Hill", .RNG.seed = 2
)

opt_naive <- run.jags(
  model = naivemodel, monitor = monitor, data = data, n.chains = 2, adapt = 5000, sample = 10000,
  burnin = 5000, method = "rjags", inits = list(inits1, inits2)
)

out_naive <- summary(opt_naive)[, c(2, 4, 5, 11)]
