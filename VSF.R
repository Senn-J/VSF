
##================================================================================##
# load Packages
##================================================================================##

library(lidR)
library(raster)
library(magrittr)
library(tidyr)
library(data.table)
library(RColorBrewer)


cellsz <- 5 # adjust cellsize for the VSF calculation



filem <- "*.las" # enter filename


# read in the .las file
lidar <- lidR::readLAS(filem)



gf_adjusted <- function(z, dz = 1, z0 = 2) #adjusted version of the lidR function gap_fraction_profile()
{
  z <- c(-0.01,z) # this makes sure all cells contain points below vegetation
  bk = seq(floor((min(z) - z0)/dz) * dz + z0, ceiling((max(z) - 
                                                         z0)/dz) * dz + z0, dz)
  if (length(bk) <= 1) 
    return(data.frame(z = numeric(0), gf = numeric(0)))
  histogram = graphics::hist(z, breaks = bk, plot = F)
  h = histogram$mids
  p = histogram$counts/sum(histogram$counts)
  p = c(p, 0)
  cs = cumsum(p)
  i = data.table::shift(cs)/cs
  i[is.na(i)] = 0
  i[is.nan(i)] = NA
  z = h
  i = i[-length(i)]
  return(data.frame(z = z[z > z0], gf = i[z > z0]))
}



# calculate metric: gap fraction profile
# this is not 2d gap fraction but 3d gap fraction as voxels
# the adjusted method adds a dummy value below the threshold to make sure values for 
#the whole profile are calculated
met <- grid_metrics3d(lidar, ~gf_adjusted(Z, dz = 1, z0 = 0.1), res = c(cellsz, 50))


#d elete redundant columns
met$Z <- met$z 
met <- met[ ,c(1, 2, 3, 5)]
#plot(met, color="lad")


# adjust table of calculated gap fraction profiles
# the four column table X-Y-Z-VALUE is transposed by Z values
# each layer z gets its own column
gaptab <- spread(met, Z, gf, drop = T)
gaptab <- data.table(gaptab) # reformat into data.table format
gaptab[is.na(gaptab)] <- 1 # assumption: all cells in gap fraction table with NA
# can be seen as 100% gap. therefore can be set to 1 (all gap)

# make sure its 50 layers 
# adding empty layers up to the table
for (i in (1 + ncol(gaptab)):52) gaptab[,as.character(i - 2.4) := 1]


# everything thats not gap is vegetation cover. vegetation cover as table
# calculating the cover as 1 - gap fraction
vegtab <- 1 - gaptab
vegtab[ ,1:2] <- gaptab[ ,1:2] # paste actual coordinates into table

#save total number of layers for later
nLay <- ncol(vegtab) - 2

#convert table to raster
gaps <- rasterFromXYZ(gaptab)

#calculate cover raster from gap raster
cov <- 1 - gaps


# now were calculating for every layer how much is intercepted in this layer or a layer below 
# the cumulative interception ratio (CIR)
CIR <- vegtab[ ,1:3]
for (i in 4:(nLay + 2)) CIR[ ,eval(names(gaptab)[i]) := vegtab[[i]] + gaptab[[i]]*CIR[[i - 1]]]
CIR[ ,open := 1]

# calculating the ratio of drips that fall to the ground without being intercepted below
# the drip contribution profile (DC)
contrib <- CIR[,1:3]
for (i in 4:(nLay + 3)) contrib[, eval(names(CIR)[i]) := -CIR[[i - 1]] + CIR[[i]]]



# make rasters, plot
interR <- rasterFromXYZ(CIR)  
contribR <- rasterFromXYZ(contrib) 
plot(interR)
plot(contribR)


# calcula the VSF
# the drip contribution profile is weighted by "weigh" the weighting factor
# and summed up for every gridcell
weigh <- c(3 - 3 * exp(-0.4*seq(0.5,nLay)),1)
vsfT <- cbind(contrib[,1:2],rowSums(t(t(contrib[,3:ncol(contrib)])*weigh)))
vsfR <- rasterFromXYZ(vsfT)

# plot result
plot(vsfR)