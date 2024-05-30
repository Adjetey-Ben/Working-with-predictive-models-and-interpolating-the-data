

#########################################################################################################################
library(rgdal)
library(tmap)
library(spTransform)
library(spatstat)  
library(maptools)  
library(raster)
library(sp)    
library(rgeos)
library(tmaptools)
library(gstat)
setwd("/home/rscc/Downloads/Data")
data<-read.table("Tmean_SSP126.txt" ,header=TRUE, na.string="-99.9", sep="",dec = ".")
df = data.frame(as.matrix(data))
coordinates(df) <- ~Lon+Lat
proj4string(df)=CRS("+init=epsg:4326") # set it to lat-long
pts = spTransform(df,CRS(proj4string(df)))
saveRDS(pts, "pts.rds", ascii=F)
p <- readRDS("pts.rds")

crswgs84=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
z=readShapePoly("/home/rscc/Downloads/North_idw_interpolation/Shape_file/GHA_adm0.shp",proj4string=crswgs84,verbose=TRUE)
saveRDS(z, "z.rds", ascii=F)

w <- readRDS("z.rds")
p@bbox <- w@bbox

png("Tmean_SSP126__Location_Stotal.png")

tm_shape(w) + tm_polygons() +
  tm_shape(p) +
  tm_dots(col="Data", palette = "RdBu", auto.palette.mapping = FALSE,
          title="Rainfall(mm)",size=0.7) +
  tm_text("Data", just="left", xmod=.5, size = 0.7) +
  tm_legend(legend.outside=TRUE)

dev.off()

#Theissen
library(spatstat)  # Used for the dirichlet tessellation function
library(maptools)  # Used for conversion from SPDF to ppp
library(raster)    # Used to clip out thiessen polygons

png("Theissen_Stotal.png")

th  <-  as(dirichlet(as.ppp(p)), "SpatialPolygons")
proj4string(th) <- proj4string(p)
th.z     <- over(th, p, fn=mean)
th.spdf  <-  SpatialPolygonsDataFrame(th, th.z)
th.clp   <- raster::intersect(w,th.spdf)

# Map the data
tm_shape(th.clp) + 
  tm_polygons(col="Data", palette="RdBu", auto.palette.mapping=FALSE,
              title="Rainfall(mm)") +
  tm_legend(legend.outside=TRUE)

dev.off()

#IDW
#library(gstat) # Use gstat's idw routine
#library(sp)    # Used for the spsample function

# Create an empty grid where n is the total number of cells
png("Tmean_SSP126_IDW_Stotal.png")
grd              <- as.data.frame(spsample(p, "regular", n=50000))
names(grd)       <- c("X", "Y")
coordinates(grd) <- c("X", "Y")
gridded(grd)     <- TRUE  # Create SpatialPixel object
fullgrid(grd)    <- TRUE  # Create SpatialGrid object

# Add P's projection information to the empty grid
proj4string(grd) <- proj4string(p)

# Interpolate the grid cells using a power value of 2 (idp=2.0)
p.idw <- gstat::idw(Data ~ 1, p, newdata=grd, idp=2.0)

# Convert to raster object then clip to Texas
r       <- raster(p.idw)
r.m     <- mask(r, w)

# Plot
tm_shape(r.m) + 
  tm_raster(n=7,palette = "OrRd", auto.palette.mapping = FALSE,
            title="Rainfall(mm)") + 
  #tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()

#Fine tuning
# Leave-one-out validation routine

png("Tmean_SSP126__Fine_tune_Stotal.png")

IDW.out <- vector(length = length(p))
for (i in 1:length(p)) {
  IDW.out[i] <- idw(Data ~ 1, p[-i,], p[i,], idp=2.0)$var1.pred
}

# Plot the differences
OP <- par(pty="s", mar=c(4,3,0,0))
plot(IDW.out ~ p$Data, asp=1, xlab="Observed", ylab="Predicted", pch=16,
     col=rgb(0,0,0,0.5))
abline(lm(IDW.out ~ p$Data), col="red", lw=2,lty=2)
abline(0,1)
par(OP)

dev.off()

# Compute RMSE
sqrt( sum((IDW.out - p$Data)^2) / length(p))

#Cross-validation
# Implementation of a jackknife technique to estimate 
# a confidence interval at each unsampled point.

# Create the interpolated surface

png("Tmean_SSP126__Confidence_I_Stotal.png")

img <- gstat::idw(Data~1, p, newdata=grd, idp=2.0)
n   <- length(p)
Zi  <- matrix(nrow = length(img$var1.pred), ncol = n)

# Remove a point then interpolate (do this n times for each point)
st <- stack()
for (i in 1:n){
  Z1 <- gstat::idw(Data~1, p[-i,], newdata=grd, idp=2.0)
  st <- addLayer(st,raster(Z1,layer=1))
  # Calculated pseudo-value Z at j
  Zi[,i] <- n * img$var1.pred - (n-1) * Z1$var1.pred
}

# Jackknife estimator of parameter Z at location j
Zj <- as.matrix(apply(Zi, 1, sum, na.rm=T) / n )

# Compute (Zi* - Zj)^2
c1 <- apply(Zi,2,'-',Zj)            # Compute the difference
c1 <- apply(c1^2, 1, sum, na.rm=T ) # Sum the square of the difference

# Compute the confidence interval
CI <- sqrt( 1/(n*(n-1)) * c1)

# Create (CI / interpolated value) raster
img.sig   <- img
img.sig$v <- CI /img$var1.pred 

# Clip the confidence raster to Texas
r <- raster(img.sig, layer="v")
r.m <- mask(r, w)

# Plot the map
tm_shape(r.m) + tm_raster(n=7,title="95% confidence interval \n(Days)") +
  #tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()

#Fit first order Polynomial
# Define the 1st order polynomial equation

png("Tmean_SSP126__Polynomial1_Stotal.png")

f.1 <- as.formula(Data ~ X + Y) 

# Add X and Y to P
p$X <- coordinates(p)[,1]
p$Y <- coordinates(p)[,2]

# Run the regression model
lm.1 <- lm( f.1, data=p)

# Use the regression model output to interpolate the surface
dat.1st <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.1, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.1st)
r.m <- mask(r, w)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=7, palette="RdBu", auto.palette.mapping=FALSE, 
            title="Rainfall(mm)") +
  #tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()

#2nd order polynomial
# Define the 2nd order polynomial equation

png("Tmean_SSP126__Polynomial2_Stotal.png")

f.2 <- as.formula(Data ~ X + Y + I(X*X)+I(Y*Y) + I(X*Y))

# Add X and Y to P
p$X <- coordinates(p)[,1]
p$Y <- coordinates(p)[,2]

# Run the regression model
lm.2 <- lm( f.2, data=p)

# Use the regression model output to interpolate the surface
dat.2nd <- SpatialGridDataFrame(grd, data.frame(var1.pred = predict(lm.2, newdata=grd))) 

# Clip the interpolated raster to Texas
r   <- raster(dat.2nd)
r.m <- mask(r, w)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=7, palette="RdBu", auto.palette.mapping=FALSE,
            title="Rainfall(mm)") +
  #tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()

############################################################################
#Kriging

#Variogram
# Define the 1st order polynomial equation

png("Tmean_SSP126__Variogram_Stotal.png")
f.1 <- as.formula(Data ~ X + Y) 

# Compute the sample variogram; note that the f.1 trend model is one of the
# parameters passed to variogram(). This tells the function to create the 
# variogram on the de-trended data.
var.smpl <- variogram(f.1, p, cloud = FALSE, cutoff=1000000, width=89900)

# Compute the variogram model by passing the nugget, sill and range values
# to fit.variogram() via the vgm() function.
dat.fit  <- fit.variogram(var.smpl, fit.ranges = FALSE, fit.sills = FALSE,
                          vgm(psill=14, model="Sph", range=790, nugget=0))

# The following plot allows us to assess the fit
plot(var.smpl, dat.fit, xlim=c(0,150))

dev.off()
####################Krigging
# Define the trend model
png("Tmean_SSP126__Krigging_Stotal.png")
f.1 <- as.formula(Data ~ X + Y) 

# Perform the krige interpolation (note the use of the variogram model
# created in the earlier step)
dat.krg <- krige( f.1, p, grd, dat.fit)

# Convert kriged surface to a raster object for clipping
r <- raster(dat.krg)
r.m <- mask(r, w)

# Plot the map
tm_shape(r.m) + 
  tm_raster(n=7, palette="OrRd", auto.palette.mapping=FALSE, 
            title="Rainfall(mm)") +
  #tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()

###Generate the variance and confidence interval maps

#The dat.krg object stores not just the interpolated values, but the variance values as well. These can be passed to the raster object for mapping as follows:

png("Tmean_SSP126__Confidence_Variance_Krig_Stotal.png")

r   <- raster(dat.krg, layer="var1.var")
r.m <- mask(r, w)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="Variance map \n(Rainfall(mm))") +tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()
#########Confidence Interval

png("Tmean_SSP126__ConfidenceI_Krig_Stotal.png")
r   <- sqrt(raster(dat.krg, layer="var1.var")) * 1.96
r.m <- mask(r, w)

tm_shape(r.m) + 
  tm_raster(n=7, palette ="Reds",
            title="95% CI map \n(Rainfall(mm))") +tm_shape(p) + tm_dots(size=0.2) +
  tm_legend(legend.outside=TRUE)

dev.off()
