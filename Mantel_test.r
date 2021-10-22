
library(dplyr)
library(ade4)
# Performing Mantel test

# Reading data from website
ozone <- read.table("https://stats.idre.ucla.edu/stat/r/faq/ozone.csv", sep=",", header=T)
head(ozone, n=10) 

# To run a Mantel test, we will need to generate two distance matrices: 
# one containing spatial distances and one containing distances between measured outcomes at the given points.  
# In the spatial distance matrix, entries for pairs of points that are close together are lower than for pairs of points that are far apart.  
# In the measured outcome matrix, entries for pairs of locations with similar outcomes are lower than for pairs of points with dissimilar outcomes.  
# We do this using the dist function.  
#The Mantel test function will require objects of this “distance” class.
station.dists <- dist(cbind(ozone$Lon, ozone$Lat))
ozone.dists <- dist(ozone$Av8top)

as.matrix(station.dists)[1:5, 1:5]
as.matrix(ozone.dists)[1:5, 1:5]

# These are the two matrices which the function will be testing for a correlation. 

mantel.rtest(station.dists, ozone.dists, nrepet = 9999)
