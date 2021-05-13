
# R code for the blogpost at
#  https://mmeredith.net/blog/2020/Model_averaging.htm

# Model averaging with wiqid::occSS output

# Need wiqid 0.3.0.9001 or later
# If necessary do
# remotes::install_github("mikemeredith/wiqid")

# Which models to include?
# ========================

library(wiqid)
data(toves)
str(toves)

# Pull out the detection data
DH <- toves[, 1:4]
mean(rowSums(DH) > 0)  # naive occupancy

# Try all 8 models
allMod <- allCombinations(response="psi", covars=c("x1", "x2", "x3"))
models <- vector("list", 8)
for(j in 1:8)
  models[[j]] <- occSS(DH, as.formula(allMod[[j]]), data=toves)
AICcs <- sapply(models, AICc)
names(AICcs) <- allMod
AICcs <- sort(AICcs)
data.frame(AICc = AICcs, delta = AICcs - min(AICcs))
#                        AICc      delta
# psi ~ x1 + x2      632.8277  0.0000000
# psi ~ x1           633.0500  0.2223035
# psi ~ x1 + x2 + x3 634.5528  1.7251486
# psi ~ x1 + x3      634.6857  1.8579871
# psi ~ 1            668.5309 35.7032110
# psi ~ x2           670.0420 37.2142704
# psi ~ x3           670.3150 37.4872842
# psi ~ x2 + x3      671.8358 39.0081545

# Models without x1 have huge deltaAICs.
# Covar x3 is also useless.

# Do model-average predictions
# ============================

# Now do the proper models
m.1 <- occSS(DH, psi ~ x1, data=toves)
m.12 <- occSS(DH, psi ~ x1 + x2, data=toves)
AICtable(AICc(m.1, m.12))
#      df    AICc Delta ModelLik ModelWt
# m.12  4 632.828 0.000    1.000   0.528
# m.1   3 633.050 0.222    0.895   0.472

# Model averaging for psi for m.1 and m.12
# --------------------------------------

# As an example, get model-averaged psi for the first 6 sites:
newdata <- toves[1:6, ]

( psi.ma <- predictAvg(list(m.1, m.12), newdata, "psi", type="response") )
#         est         SE     lowCI     uppCI
# 1 0.8515831 0.13036260 0.4318010 0.9774376
# 2 0.2937862 0.07341358 0.1721334 0.4542412
# 3 0.6673564 0.22612654 0.2141285 0.9365958
# 4 0.2748484 0.06956039 0.1605406 0.4289563
# 5 0.9134818 0.07464093 0.6238355 0.9853414
# 6 0.3079147 0.07877618 0.1773493 0.4786724

# Get individual model predictions for comparison
( psi.1 <- predict(m.1, newdata, parameter="psi", type="response") )
( psi.12 <- predict(m.12, newdata, parameter="psi", type="response") )


# Plot point estimates and CIs
plot(1:6, psi.ma[,'est'], xlab="Site number", ylab="psi", pch=16, cex=1.5, 
    las=1, ylim=0:1, xlim=c(0.5, 6.5))
arrows(1:6, psi.ma[,'lowCI'], 1:6, psi.ma[,'uppCI'], angle=90, length=0.03, code=3, lwd=2)
# Add values from psi.1 and psi.12
points(1:6 - 0.2, psi.1[,1], col='red')
arrows(1:6 - 0.2, psi.1[,3], 1:6 - 0.2, psi.1[,4],
    angle=90, length=0.03, code=3, col='red')
points(1:6 + 0.2, psi.12[,1], pch=2, col='blue')
arrows(1:6 + 0.2, psi.12[,3], 1:6 + 0.2, psi.12[,4],
    angle=90, length=0.03, code=3, col='blue')

# Model-averaged CIs are mostly wider than from the underlying models, especially
#   when the point estimates from the models are very different, eg. site 3.
# That's because it includes uncertainty about which model is best.

# Plotting relationships
# ----------------------

range(toves$x1)
# [1] 0.07 9.99
x1 <- seq(0, 10, , 100)
range(toves$x2)
# [1]  0.46 10.26
x2 <- seq(0, 10, , 100)
newdata <- expand.grid(x1, x2)
colnames(newdata) <- c("x1", "x2")
head(newdata)
#          x1 x2
# 1 0.0000000  0
# 2 0.1010101  0
# 3 0.2020202  0
# 4 0.3030303  0
# 5 0.4040404  0
# 6 0.5050505  0

# Get model averaged predictions
psi.ma <- predictAvg(list(m.1, m.12), newdata, "psi", type="response")
mat.ma <- matrix(psi.ma[, 'est'], 100, 100)

# Get predictions from the individual models
psi.1 <- predict(m.1, newdata, parameter="psi", type="response")
mat.1 <- matrix(psi.1[, 'est'], 100, 100)
psi.12 <- predict(m.12, newdata, parameter="psi", type="response")
mat.12 <- matrix(psi.12[, 'est'], 100, 100)

# Do the image plots
x11(width=7, height=2.73)
par(mfrow = c(1,3))
image(x1, x2, mat.1, asp=1, main="psi ~ x1")
box()
image(x1, x2, mat.ma, asp=1, main="Model average psi")
box()
image(x1, x2, mat.12, asp=1, main="psi ~ x1 + x2")
box()
