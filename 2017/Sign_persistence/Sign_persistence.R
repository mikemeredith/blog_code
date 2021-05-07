
# Code to create the plots displayed in this blog post:
#  https://mmeredith.net/blog/2017/Sign_persistence.htm

source("[function]_MIKE_persistence.R")

dung <- read.csv("Nakai_dung_disapp.csv")

summary(dung)
dim(dung)
plot(STATE ~ DAYS, data=dung)

# Dung piles were classified as 'visible' (STATE = 1) or 'disappeared' (STATE = 0).
# DAYS is the number of days *prior* to the density survey that the dung was deposited.

# Code for the first plot, showing the data
# -----------------------------------------
# Since DAYS = time prior to Day 0, we'll plot it from -500 to 0 instead of
# from 0 to 500 days. We'll also spread out the points (with jittering) and
# improve the labelling.

xlim <- range(-dung$DAYS, 0)
plot(jitter(-dung$DAYS), jitter(dung$STATE), xlim=xlim, pch=18, cex=0.3,
   xlab="Time of deposition\n(days before density survey)", mgp=c(2.3,0.3,0),
   tcl=NA,
   ylab="State of dung\nat time of density survey", yaxt='n')
abline(h=0:1, lty=2)
axis(2, at=0:1, label=c("disappeared", "visible"), tcl=NA)
abline(v=0, col='blue', lwd=2)
mtext("Dung density surveys centered on Day 0", adj=1, col='blue')

# Second plot, showing the estimate
# ---------------------------------

# Do the necessary calculations

# The extimate
( EST <- MIKE.persistence(dung) )

# The decay curve
dung.glm <- glm(STATE ~ DAYS, data=dung, family=binomial(link = "logit"))
x <- seq(0, max(dung$DAYS), length=100)
pred <- predict(dung.glm, newdata=list(DAYS=x), type='response')

# Redraw the plot _without_ the points (type = 'n'):
plot(jitter(-dung$DAYS), jitter(dung$STATE), xlim=xlim, pch=18, cex=0.3,
   xlab="Time of deposition\n(days before density survey)", mgp=c(2.3,0.3,0),
   tcl=NA, type='n',
   ylab="State of dung\nat time of density survey", yaxt='n')
axis(2, at=0:1, label=c("disappeared", "visible"), tcl=NA)
abline(v=0, col='blue', lwd=2)
mtext("Dung density surveys centered on day 0", adj=1, col='blue')

# Draw the shaded polygons next, or they will hide the points:
est <- round(EST[1])
x1 <- seq(0, est, length=50)
pred1 <- predict(dung.glm, newdata=list(DAYS=x1), type='response')
x1 <- c(est, x1) ; pred1 <- c(1, pred1)
polygon(-x1, pred1, col='lightpink', border=NA)
x2 <- seq(est, max(dung$DAYS), length=50)
pred2 <- predict(dung.glm, newdata=list(DAYS=x2), type='response')
x2 <- c(x2,est) ; pred2 <- c(pred2, 0)
polygon(-x2, pred2, col='lightpink', border=NA)

# Now add the red lines and explanations:
lines(-x, pred, col='red', lwd=2)
abline(v=-est, lty=3, col='red')
text(-est/2, 0.3, paste("survival =\n", round(est), "days"), col='red')
text(-200, 0.8, "Some dung <150 days old disappears", col='red', pos=2, cex=0.9)
arrows(-200, 0.8, -140, 0.8, length=0.1, col='red')
text(-200, 0.3, "Some dung >150 days old is still visible", col='red', pos=2, cex=0.9)
arrows(-200, 0.3, -160, 0.2, length=0.1, col='red')
text(-180, 0.5, "Curve shows probability of surviving", pos=2, col='red')

# Finally, add the points and the horizontal dotted lines:
points(jitter(-dung$DAYS), jitter(dung$STATE), xlim=xlim, pch=18, cex=0.3)
abline(h=0:1, lty=3)

# Print it out and hang it on your office wall!
