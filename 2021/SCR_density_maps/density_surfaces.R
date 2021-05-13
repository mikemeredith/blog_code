
# R code for blog post at
#  https://mmeredith.net/blog/2021/SCR_density_maps.htm

# Density surfaces
# ================

# Use example in makeJAGSmask package
# If necessary
# remotes::install_github("mikemeredith/makeJAGSmask")

library(makeJAGSmask)
data(simSCR)
str(simSCR, 2)

# Get the bits we need
# --------------------

sigma <- simSCR$sims.list$sigma
S <- simSCR$sims.list$S
w <- simSCR$sims.list$w
S[w == 0] <- NA  # replace coordinates with NA for phantoms

habMat <- simSCR$JAGSmask$habMat
trapMat <- simSCR$JAGSmask$trapMat
upperLimit <- simSCR$JAGSmask$upperLimit

niter <- dim(w)[1]

# Plot the MCMC chains for the activity centres
# ---------------------------------------------
windows(4.74, 5)
image(x=(1:65)+0.5, y=(1:59)+0.5, z=habMat, asp=1, xlab="", ylab="", col=c("grey", "white"))
for(i in 1:5)
  points(S[,i,], cex=0.1, col=i+1)
for(i in 6:30)
  points(S[,i,], cex=0.1, col=1)
points(trapMat, pch=3, col='red')

# Posterior probability density
# =============================

# A useful function to "tot up" the number of items in each pixel
# index = a 3D array, iterations x species x 2 with indices into the matrix 'mat'
#   non-integers are truncated, NAs are ignored
# the specified elements in the matrix will be augmented by 1 for each occurrence in 'index'
totup <- function(index, mat) {
  index <- matrix(index, length(index)/2, 2)
  for(i in 1:nrow(index))
    mat[index[i,1], index[i,2]] <- mat[index[i,1], index[i,2]] + 1
  return(mat)
}


# posterior probability density map for the AC of 1 animal
ACmat <- array(0, dim=dim(habMat))
ACmat1 <- totup(S[,2,], ACmat)
ACmat <- ACmat1 / niter
max(ACmat)
ACmat[habMat == 0] <- NA
ACcol <- heat.colors(12, rev=TRUE)
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ACmat, asp=1, xlab="", ylab="", col=ACcol)
points(trapMat, pch=3, cex=0.5)
# Zoom in
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ACmat, asp=1, xlab="", ylab="", col=ACcol,
    xlim=c(3,13), ylim=c(20,30))
points(trapMat, pch=3, cex=0.5)

# Activity centre density
# -----------------------

# Now aggregate over all animals
ACmat <- array(0, dim=dim(habMat))
ACmat1 <- totup(S, ACmat)
ACmat <- ACmat1 / niter
ACmat[habMat == 0] <- NA
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ACmat, asp=1, xlab="", ylab="", col=ACcol)
points(trapMat, pch=3, cex=0.5)

# Animal density
# ==============

set.seed(42)
# Array for the posterior predictive animal locations
AL <- array(NA, dim=dim(S))
nind <- dim(S)[2]
# Work through all iterations and all animals
for(i in 1:niter)
  for(j in 1:nind)
    if(w[i,j] == 1) {
      hab <- 0
      while(hab == 0) {
        AL[i,j, ] <- S[i,j, ] + rnorm(2, 0, sigma[i])
        # habitat check
        if(all(AL[i,j,] >= 1) && all(AL[i,j,] < upperLimit))
          hab <- habMat[trunc(AL[i,j,1]),trunc(AL[i,j,2])]
      }
    }

# posterior predictive animal density map for 1 animal
ALmat <- array(0, dim=dim(habMat))
ALmat <- totup(AL[,2,], ALmat)
ALmat <- ALmat / niter
ALmat[habMat == 0] <- NA
ALcol <- terrain.colors(12, rev=TRUE)
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ALmat, asp=1, xlab="", ylab="", col=ALcol)
points(trapMat, pch=3, cex=0.5)
# Zoom in
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ALmat, asp=1, xlab="", ylab="", col=ALcol,
    xlim=c(1,16), ylim=c(17,32))
points(trapMat, pch=3)

# Now aggregate over all animals
ALmat <- array(0, dim=dim(habMat))
ALmat <- totup(AL, ALmat)
ALmat <- ALmat / niter
ALmat[habMat == 0] <- NA
image(x=(1:65)+0.5, y=(1:59)+0.5, z=ALmat, asp=1, xlab="", ylab="", col=ALcol)
points(trapMat, pch=3, cex=0.5)

Dmat <- ALmat * 100 # convert to density per 100 km2

# Do more simulated animal locations
set.seed(42)
# Array for the posterior predictive animal locations
AL <- array(NA, dim=dim(S))
ALmat <- array(0, dim=dim(habMat))
nind <- dim(S)[2]
# Work through all iterations and all animals FIVE TIMES
for(q in 1:5) {
  for(i in 1:niter)
    for(j in 1:nind)
      if(w[i,j] == 1) {
        hab <- 0
        while(hab == 0) {
          AL[i,j, ] <- S[i,j, ] + rnorm(2, 0, sigma[i])
          # habitat check
          if(all(AL[i,j,] >= 1) && all(AL[i,j,] < upperLimit))
            hab <- habMat[trunc(AL[i,j,1]),trunc(AL[i,j,2])]
        }
      }

  ALmat <- totup(AL, ALmat)
}
ALmat <- ALmat / niter / 5
ALmat[habMat == 0] <- NA

Dmat <- ALmat * 100 # convert to density per 100 km2

# Converting to a raster
# ======================
library(raster)

# Starting from a data frame of coordinates
# -----------------------------------------
gridPix <- cbind(c(row(Dmat))+0.5, c(col(Dmat))+0.5) # coordinates on pixel scale
head(gridPix)
gridM <- convertOutput(gridPix, simSCR$JAGSmask) # convert to real scale
head(gridM)
xyz <- cbind(gridM, c(Dmat))
densityR <- rasterFromXYZ(xyz)
plot(densityR)
points(simSCR$traps, pch=3, cex=0.5)
points(simSCR$trueACs, pch=16)

# Starting from a habitat raster
# ------------------------------
plot(simSCR$patchR)
points(simSCR$traps, pch=3, col='red')

densityR <- simSCR$patchR  # copy the original habitat raster
tmp <- t(Dmat)             # transpose the matrix...
tmp <- tmp[nrow(tmp):1,]   # ...then flip vertically
values(densityR) <- tmp
st <- stack(simSCR$patchR, densityR)
names(st) <- c("edge_distance", "density")
plot(st)
