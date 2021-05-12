
# R code for the blog post at
#  https://mmeredith.net/blog/2020/storing_MCMC.htm

# If necessary:
# install.packages("mcmcOutput")
library(mcmcOutput)
data(mcmcListExample)
str(mcmcListExample)

# convert to class mcmcOutput
mco <- mcmcOutput(mcmcListExample)
str(mco)
print(mco)
summary(mco)
View(summary(mco))

mcos <- sumryList(mco)
mcos$mean$p
mcos$MCEpc$psi

# Extract with "$"
p <- mco$p
str(p)
p[1:5,,]  # Elements of p not defined in the model are filled with NAs

# "[" with one index
head(mco[4:5])
head(mco[c("z[35]", "z[39]")])

# "[" with two indices
mco[1:5, "psi"]  # First 5 values for psi (chain #1)

# "[" with three indices
mco[1:5, 2, "psi"] # First 5 values for psi in chain #2

# Plotting methods/functions
diagPlot(mco)
diagPlot(mcmcListExample) # using the original mcmc.list object
plot(mco)
plot(mco, compVal=0) # show proportion above/below 0.

 

