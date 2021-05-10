
# Code for blog post at
#  https://mmeredith.net/blog/2019/makeJAGSmask_enhancements.htm

# Youmay need to do this
# install.packages(c("abind", "plotrix", "raster", "secr", "rgeos"))
# install.packages("remotes")
# remotes::install_github("mikemeredith/makeJAGSmask")
packageVersion("makeJAGSmask")  # 0.1.1 or later

# Using 'core'
# ------------

library(makeJAGSmask)
library(sp)
# Get the example data and plot them:
data(simSCR)
str(simSCR, 1)
attach(simSCR)
plot(JAGSmask)

# Convert AC locations to the original units, remove phantoms:
AC <- convertOutput(sims.list$S, JAGSmask)
AC[sims.list$w == 0] <- NA
plotACs(which=1:30, AC=AC, traps=traps, Y=Y, hab=patchDF, show.labels=FALSE)

# There is a river flowing through the study area, and the portion north of the river lies in the municipality of Arkadia.
polygon(Arkadia, border='red')

coremask <- addCore(JAGSmask, type="poly", poly=Arkadia)
plot(coremask)

N <- getACinCore(sims.list$S, sims.list$w, coremask)
plot(table(N)) # posterior probability distribution
mean(N)        # posterior mean
quantile(N, probs=c(0.025, 0.975))  # 95% credible interval

# Using rasters
# -------------
library(raster)
plot(patchRS)

edge <- patchRS$edge
patchRS$edgeS <- (edge - mean(values(edge), na.rm=TRUE)) / sd(values(edge), na.rm=TRUE)
patchRS

# Now we create our JAGSmask object which will include the raster information as matrices:
mymask <- convertRaster(patchRS, traps)
str(mymask)

M <- 25  # size of the augmented data set
nTraps <- nrow(traps)
yAug <- rbind(Y, matrix(0, M - nrow(Y), nTraps))

# organise data for JAGS
dataList <- c(mymask, list(y = yAug, M = M, nTraps = nTraps, zeros = rep(0, M), nOcc = 90))
str(dataList)

# Initial values
inits <- function() list(w = rep(1, M), S = matrix(15, M, 2))

wanted <- c("N", "omega", "p0", "beta", "sigma", "S", "w")
library(jagsUI)
result <- jags(dataList, inits, wanted, "patch_covar.jags", DIC=FALSE,
  n.chains=3, n.adapt=1000, n.iter=5500, n.burnin=500, n.thin=10, parallel=TRUE)
result
