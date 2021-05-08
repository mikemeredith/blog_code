
# Code for the blog post at
#  https://mmeredith.net/blog/2018/Colours.htm

# Colour-blind-friendly palette
# -----------------------------

palette("default")  # Use the default palette
barplot(rep(1,10), yaxt="n", col=1:10, names.arg=1:10)

library(dichromat)
palette(dichromat(palette()))
barplot(rep(1,10), yaxt="n", col=1:10, names.arg=1:10)

# Here is the palette of colour-blind friendly colours proposed by Okabe and Ito:
.cbColors <- c( black = "#000000", # black
        vermilion = "#D55E00",     # red
        green = "#009E73",         # green3
        blue = "#0072B2",          # blue
        sky = "#56B4E9",           # cyan
        purple = "#CC79A7",        # magenta
        yellow = "#F0E442",        # yellow
        orange = "#E69F00")        # (grey)
grDevices::palette(.cbColors)

palette(.cbColors)  # Use CB-friendly palette
barplot(rep(1,10), yaxt="n", col=1:10, names.arg=1:10)

library(dichromat)
palette(dichromat(palette()))
barplot(rep(1,10), yaxt="n", col=1:10, names.arg=1:10)


# Point size and shape, line width and style
# ------------------------------------------

# Use the first 30 rows of the 'trees' data set for the scatter plots
tr <- trees[1:30, ]

par(mfrow=c(1,2))
palette("default")
with(tr, scatter.smooth(Height, Volume, lpars=list(col=2, lwd=1), bty='l',
col=1:4,
  main="Small, narrow"))
abline(lm(Volume ~ Height, data=trees), lwd=1, col=3)
legend('topleft', c("loess", "regression"), lwd=1, col=2:3, bty='n')

palette(.cbColors)
with(tr, scatter.smooth(Height, Volume, lpars=list(col=2, lwd=4), bty='l',
  col=1:3, pch=15:17, main="Big,  wide"))
abline(lm(Volume ~ Height, data=trees), lwd=4, col=3, lty=3)
legend('topleft', c("loess", "regression"), lwd=4, col=2:3, bty='n', lty=c(1,3))


palette(.cbColors)
par(mar=c(5,4,4,6)+.1)
with(tr, scatter.smooth(Height, Volume, lpars=list(col=2, lwd=4), bty='l', las=1,
  main="Labelling on the plot"))
abline(lm(Volume ~ Height, data=trees), lwd=4, col=3, lty=3)
text(86, 32, "loess", pos=4, xpd=TRUE)
text(86, 48, "regression", pos=4, xpd=TRUE)

# ... or this
par(mar=c(5,4,4,2)+.1)
with(tr, scatter.smooth(Height, Volume, lpars=list(col=2, lwd=4), bty='l', las=1,
  main="Labelling on the plot"))
abline(lm(Volume ~ Height, data=trees), lwd=4, col=3, lty=3)
text(78, 20, 'loess', pos=1) 
arrows(78, 19.58, 76.528, 26.658, length=0.1)
text(67, 34, 'regression', pos=3) 
arrows(67, 34, 72.87, 25.73, length=0.1) 
