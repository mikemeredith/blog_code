
# Code for the blog post at
#  https://mmeredith.net/blog/2017/Arranging_arrays.htm

colNames <- as.vector(t(outer(1:2, 1:4, paste, sep=".")))
siteNames <- c("A", "B", "C")
dat1 <- matrix(rpois(24, 1.5), 3)
dimnames(dat1) <- list(site = siteNames, occasion = colNames)
dat1
#     occasion
# site 1.1 1.2 1.3 1.4 2.1 2.2 2.3 2.4
#    A   3   1   3   3   2   2   1   3
#    B   0   1   0   0   2   2   1   5
#    C   3   0   0   0   2   2   3   4

mat1 <- outer(siteNames, colNames, paste0)
dimnames(mat1) <- list(site = siteNames, occasion = colNames)
mat1
#     occasion
# site 1.1    1.2    1.3    1.4    2.1    2.2    2.3    2.4
#    A "A1.1" "A1.2" "A1.3" "A1.4" "A2.1" "A2.2" "A2.3" "A2.4"
#    B "B1.1" "B1.2" "B1.3" "B1.4" "B2.1" "B2.2" "B2.3" "B2.4"
#    C "C1.1" "C1.2" "C1.3" "C1.4" "C2.1" "C2.2" "C2.3" "C2.4"

# A comment on dimnames
# ---------------------

mat1[, 1:4]
#     occasion
# site 1.1    1.2    1.3    1.4
#    A "A1.1" "A1.2" "A1.3" "A1.4"
#    B "B1.1" "B1.2" "B1.3" "B1.4"
#    C "C1.1" "C1.2" "C1.3" "C1.4"

mat1[2, ]
#    1.1    1.2    1.3    1.4    2.1    2.2    2.3    2.4
# "B1.1" "B1.2" "B1.3" "B1.4" "B2.1" "B2.2" "B2.3" "B2.4"

dim(mat1[2, ])    # no dimensions, it's not a matrix
# NULL

dim(mat1[2, , drop = FALSE])    # now okay
# [1] 1 8

mat1[2, , drop = FALSE]    # now okay
#     occasion
# site 1.1    1.2    1.3    1.4    2.1    2.2    2.3    2.4
#    B "B1.1" "B1.2" "B1.3" "B1.4" "B2.1" "B2.2" "B2.3" "B2.4"

# Convert our simple matrix to a 3-D array
# ----------------------------------------

str(mat1)
# chr [1:3, 1:8] "A1.1" "B1.1" "C1.1" "A1.2" "B1.2" "C1.2" "A1.3" "B1.3" ...
# - attr(*, "dimnames")=List of 2
#  ..$ site    : chr [1:3] "A" "B" "C"
#  ..$ occasion: chr [1:8] "1.1" "1.2" "1.3" "1.4" ...

# Make an array with sites (rows) x visits (columns) x years (pages)

( arr1 <- array(mat1, c(3, 4, 2)) )
# , , 1

#      [,1]   [,2]   [,3]   [,4]
# [1,] "A1.1" "A1.2" "A1.3" "A1.4"
# [2,] "B1.1" "B1.2" "B1.3" "B1.4"
# [3,] "C1.1" "C1.2" "C1.3" "C1.4"
#
# , , 2
#
#      [,1]   [,2]   [,3]   [,4]
# [1,] "A2.1" "A2.2" "A2.3" "A2.4"
# [2,] "B2.1" "B2.2" "B2.3" "B2.4"
# [3,] "C2.1" "C2.2" "C2.3" "C2.4"

dimnames(arr1) <- list(site = siteNames, visit=1:4, year=1:2)
arr1
# , , year = 1
#
#     visit
# site 1      2      3      4
#    A "A1.1" "A1.2" "A1.3" "A1.4"
#    B "B1.1" "B1.2" "B1.3" "B1.4"
#    C "C1.1" "C1.2" "C1.3" "C1.4"
#
# , , year = 2
#
#     visit
# site 1      2      3      4
#    A "A2.1" "A2.2" "A2.3" "A2.4"
#    B "B2.1" "B2.2" "B2.3" "B2.4"
#    C "C2.1" "C2.2" "C2.3" "C2.4"

arr1[3,,] # site 3
#      year
# visit 1      2
#     1 "C1.1" "C2.1"
#     2 "C1.2" "C2.2"
#     3 "C1.3" "C2.3"
#     4 "C1.4" "C2.4"

arr1[,2,] # visit 2
#     year
# site 1      2
#    A "A1.2" "A2.2"
#    B "B1.2" "B2.2"
#    C "C1.2" "C2.2"
   
arr1[,,1] # year 1
#     visit
# site 1      2      3      4
#    A "A1.1" "A1.2" "A1.3" "A1.4"
#    B "B1.1" "B1.2" "B1.3" "B1.4"
#    C "C1.1" "C1.2" "C1.3" "C1.4"

# A more complicated example
# --------------------------

rowNames <- as.vector(t(outer(c("a", "b"), siteNames, paste0)))
mat2 <- outer(rowNames, colNames, paste0)
dimnames(mat2) <- list(rowNames, colNames)
mat2
#    1.1     1.2     1.3     1.4     2.1     2.2     2.3     2.4
# aA "aA1.1" "aA1.2" "aA1.3" "aA1.4" "aA2.1" "aA2.2" "aA2.3" "aA2.4"
# aB "aB1.1" "aB1.2" "aB1.3" "aB1.4" "aB2.1" "aB2.2" "aB2.3" "aB2.4"
# aC "aC1.1" "aC1.2" "aC1.3" "aC1.4" "aC2.1" "aC2.2" "aC2.3" "aC2.4"
# bA "bA1.1" "bA1.2" "bA1.3" "bA1.4" "bA2.1" "bA2.2" "bA2.3" "bA2.4"
# bB "bB1.1" "bB1.2" "bB1.3" "bB1.4" "bB2.1" "bB2.2" "bB2.3" "bB2.4"
# bC "bC1.1" "bC1.2" "bC1.3" "bC1.4" "bC2.1" "bC2.2" "bC2.3" "bC2.4"

as.vector(mat2)
 [1] "aA1.1" "aB1.1" "aC1.1" "bA1.1" "bB1.1" "bC1.1" "aA1.2" "aB1.2" "aC1.2" "bA1.2" "bB1.2"
[12] "bC1.2" "aA1.3" "aB1.3" "aC1.3" "bA1.3" "bB1.3" "bC1.3" "aA1.4" "aB1.4" "aC1.4" "bA1.4"
[23] "bB1.4" "bC1.4" "aA2.1" "aB2.1" "aC2.1" "bA2.1" "bB2.1" "bC2.1" "aA2.2" "aB2.2" "aC2.2"
[34] "bA2.2" "bB2.2" "bC2.2" "aA2.3" "aB2.3" "aC2.3" "bA2.3" "bB2.3" "bC2.3" "aA2.4" "aB2.4"
[45] "aC2.4" "bA2.4" "bB2.4" "bC2.4"

arr2 <- array(mat2, c(3, 2, 4, 2))
dimnames(arr2) <- list(site = siteNames,
                     species = c("a", "b"),
                     visit = paste0("v", 1:4),
                     year = paste0("y", 1:2))
arr2
# , , visit = v1, year = y1

#     species
# site a       b
#    A "aA1.1" "bA1.1"
#    B "aB1.1" "bB1.1"
#    C "aC1.1" "bC1.1"
#
# , , visit = v2, year = y1
#
#     species
# site a       b
#    A "aA1.2" "bA1.2"
#    B "aB1.2" "bB1.2"
#    C "aC1.2" "bC1.2"
#
# , , visit = v3, year = y1
#
#     species
# site a       b
#    A "aA1.3" "bA1.3"
#    B "aB1.3" "bB1.3"
#    C "aC1.3" "bC1.3"
#
# , , visit = v4, year = y1
#
#     species
# site a       b
#    A "aA1.4" "bA1.4"
#    B "aB1.4" "bB1.4"
#    C "aC1.4" "bC1.4"
#
# , , visit = v1, year = y2
#
#     species
# site a       b
#    A "aA2.1" "bA2.1"
#    B "aB2.1" "bB2.1"
#    C "aC2.1" "bC2.1"
#
# , , visit = v2, year = y2
#
#     species
# site a       b
#    A "aA2.2" "bA2.2"
#    B "aB2.2" "bB2.2"
#    C "aC2.2" "bC2.2"
#
# , , visit = v3, year = y2
#
#     species
# site a       b
#    A "aA2.3" "bA2.3"
#    B "aB2.3" "bB2.3"
#    C "aC2.3" "bC2.3"
#
# , , visit = v4, year = y2
#
#     species
# site a       b
#   A "aA2.4" "bA2.4"
#    B "aB2.4" "bB2.4"
#    C "aC2.4" "bC2.4"

# use the aperm() function to permute the dimensions.

( arr3 <- aperm(arr2, c(1,3,4,2)) )
# , , year = y1, species = a
#
#     visit
# site v1      v2      v3      v4
#    A "aA1.1" "aA1.2" "aA1.3" "aA1.4"
#    B "aB1.1" "aB1.2" "aB1.3" "aB1.4"
#    C "aC1.1" "aC1.2" "aC1.3" "aC1.4"
#
# , , year = y2, species = a
#
#     visit
# site v1      v2      v3      v4
#    A "aA2.1" "aA2.2" "aA2.3" "aA2.4"
#    B "aB2.1" "aB2.2" "aB2.3" "aB2.4"
#    C "aC2.1" "aC2.2" "aC2.3" "aC2.4"
#
# , , year = y1, species = b
#
#     visit
# site v1      v2      v3      v4
#    A "bA1.1" "bA1.2" "bA1.3" "bA1.4"
#    B "bB1.1" "bB1.2" "bB1.3" "bB1.4"
#    C "bC1.1" "bC1.2" "bC1.3" "bC1.4"
#
# , , year = y2, species = b
#
#     visit
# site v1      v2      v3      v4
#    A "bA2.1" "bA2.2" "bA2.3" "bA2.4"
#    B "bB2.1" "bB2.2" "bB2.3" "bB2.4"
#    C "bC2.1" "bC2.2" "bC2.3" "bC2.4"
