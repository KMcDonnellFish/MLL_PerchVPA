rm(list = ls())
load("data/dat.Rdata")

## Setup ####
ages <- unique(dat$age)
years <- unique(dat$year)
wbar <- colMeans(dat$wa)
ca <- dat$caa # catch at age for all gears

# NAs = 0
ca[which(is.na(ca))] <- 0

# Inputs
M <- 0.2
S <- exp(-M)
Ut <- 0.5
ny <- length(years)
na <- length(ages)
N <- matrix(0, nrow = ny, ncol = na) # numbers at age matrix

## VPA ####
# STEP 1:
# initialize
N[, na] <- ca[, na] / Ut
N[!is.finite(N)] <- 0 # deal with divide by zero

# STEP 2:
# run backward VPA equation using Pope's approximation
for (t in (ny - 1):1) {
  for (a in (na - 1):1) {
    N[t, a] <- N[t + 1, a + 1] * exp(M) + ca[t, a] * exp(M / 2)
  }
}

U_final <- ca / N
# gill net vulnerability

# # STEP3:
# # reinitialize last year for U and N
# UbarT <- colMeans(U_final[1:(ny-1),]) # set for new Ubar for last year
# Vbar <- UbarT/UbarT[na-1] # vulnerability
# Vbar[na] <- 1
# N[ny,] <- ca[ny,] / (Vbar * Ut) # calculate N for last year
# U_final <- ca / N


## Plotting ####
# numbers at age
col_v <- rep(c("black", "#d62728", "#2ca02c", "#1f77b4", "#00D1D1", "#9467bd"), 3)
b_labels <- ages
plot(1:ny, seq(0, max(N, na.rm = TRUE) * 1.1, length.out = ny),
  type = "n", xlim = c(1, ny),
  xlab = "Years", ylab = "Numbers at age", xaxt = "n", yaxt = "n"
)
for (j in 1:na) {
  lines(1:ny, N[, j],
    type = "b", pch = "", col = col_v[j]
  )
  text(1:ny,
    N[, j],
    label = b_labels[j],
    cex = 0.9, col = col_v[j]
  )
}
axis(2)
axis(1, at = 1:ny, labels = years)

plot(rowSums(N), type = "l")
plot(rowSums(ca), type = "l")

# cumulative abundance at age
y <- N[, 1]
x <- years
plot(x, y,
  type = "b", pch = "",
  ylim = c(0, max(rowSums(N)) * 1.1), col = 1,
  xaxt = "n", xlab = "Years", ylab = "Cumulative Numbers at age"
)
axis(1, at = x, labels = years)
text(x, y,
  label = b_labels[1],
  cex = 0.9, col = col_v[1]
)
for (a in 2:na) {
  y <- y + N[, a]
  lines(x, y, type = "b", pch = "", col = col_v[a])
  text(x, y,
    label = b_labels[a],
    cex = 0.9, col = col_v[a]
  )
}

