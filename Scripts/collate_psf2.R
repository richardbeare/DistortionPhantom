library(tidyverse)
library(oro.nifti)

# modified version that doesn't average PSF images
# Means we can fit the gaussians individually

# sag_fast - 5x1.6x1.6
# sag - 5x1.5x1.5
# cor_fast - 1.6x5x1.6
# cor - 1.5x5.1.5
# axi_fast - 1.6x1.6x5
# axi - 1.5x1.5x5
# focus on the data from calibre

extractPSFIm <- function(psf, ScanType) {
  d <- dim(psf)
  
  md <- d/2 + 1
  bb <- 25
  bbpo <- bb + 1
  cc <- 2*bb + 1
  begin <- md - bb
  theend <- md + bb
  
  middle <- psf[begin[1]:theend[1], begin[2]:theend[2], begin[3]:theend[3]]
  coords <- which(!is.na(middle), arr.ind = TRUE)
  coords <- as_tibble(coords)
  coords <- mutate(coords, brightness=as.vector(middle))
  coords <- mutate(coords,
                   brightness=brightness/max(brightness),
                   dim1 = 0.5*(dim1 - bb - 1),
                   dim2 = 0.5*(dim2 - bb - 1),
                   dim3 = 0.5*(dim3 - bb - 1))
  if (str_detect(ScanType, "AXI")) {
    coords <- rename(coords, inplaneA=dim1, inplaneB=dim2, thruplane=dim3)
  }
  if (str_detect(ScanType, "COR")) {
    coords <- rename(coords, inplaneB=dim1, thruplane=dim2, inplaneA=dim3)
  }
  if (str_detect(ScanType, "SAG")) {
    coords <- rename(coords, thruplane=dim1, inplaneB=dim2, inplaneA=dim3)
  }
  
  return(coords)
}

extractPSF <- function(psf, ScanType) {
  d <- dim(psf)
  
  md <- d/2 + 1
  bb <- 35
  bbpo <- bb + 1
  cc <- 2*bb + 1
  begin <- md - bb
  theend <- md + bb
  
  middle <- psf[begin[1]:theend[1], begin[2]:theend[2], begin[3]:theend[3]]
  if (str_detect(ScanType, "AXI")) {
    dfA <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,bbpo,], direction="through plane", type="observed")
    dfA <- mutate(dfA, value = value/max(value))
    dfB <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,,bbpo], direction="inplaneA", type="observed")
    dfB <- mutate(dfB, value = value/max(value))
    dfC <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[,bbpo,bbpo], direction="inplaneB", type="observed")
    dfC <- mutate(dfC, value = value/max(value))
    
  }
  if (str_detect(ScanType, "COR")) {
    dfA <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,bbpo,], direction="inplaneA", type="observed")
    dfA <- mutate(dfA, value = value/max(value))
    dfB <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,,bbpo], direction="through plane", type="observed")
    dfB <- mutate(dfB, value = value/max(value))
    dfC <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[,bbpo,bbpo], direction="inplaneB", type="observed")
    dfC <- mutate(dfC, value = value/max(value))
    
  }
  if (str_detect(ScanType, "SAG")) {
    dfA <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,bbpo,], direction="inplaneB", type="observed")
    dfA <- mutate(dfA, value = value/max(value))
    dfB <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,,bbpo], direction="inplaneA", type="observed")
    dfB <- mutate(dfB, value = value/max(value))
    dfC <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[,bbpo,bbpo], direction="through plane", type="observed")
    dfC <- mutate(dfC, value = value/max(value))
    
  }
  if (str_detect(ScanType, "Fast")) {
    sig_inplane <- (1.2 * 1.6)/(2*sqrt(2*log(2)))
    sig_thruplane <- 5/(2*sqrt(2*log(2)))
  } else {
    sig_inplane <- (1.2 * 1.5)/(2*sqrt(2*log(2)))
    sig_thruplane <- 5/(2*sqrt(2*log(2)))
  }
  
  gaussian <- function(X, sigma, mu=0) {
    1.0/(sigma*sqrt(2*pi)) * exp(-0.5 * (X - mu)^2/sigma^2)
  }
  
  dftheory <- tibble( xcoord = 0.5*((1:cc) - bb))
  dftheoryA <- mutate(dftheory, value = gaussian(xcoord, sig_inplane), value=value/max(value), direction="inplaneA", type="recommended")
  dftheoryB <- mutate(dftheory, value = gaussian(xcoord, sig_inplane), value=value/max(value), direction="inplaneB", type="recommended")
  dftheoryC <- mutate(dftheory, value = gaussian(xcoord, sig_thruplane), value=value/max(value), direction="through plane", type="recommended")
  
  df <- bind_rows(dftheoryA, dftheoryB, dftheoryC, dfA, dfB, dfC)
  df <- mutate(df, ScanType = ScanType)
  return(df)
}
# 
fl <- "calibre_psf2.Rda"
# now do an iterative approach - first pass with the psf from literature applied to the target image.
# after that, replace with the newly estimated values
#fl <- "calibre_psf2_B.Rda"



# this one is unwarped first - seems to give crazy results
#fl <- "calibre_psf2_uw.Rda"
if (file.exists(fl)) {
  load(fl)
} else {
  calibre <- tibble(filename=readLines("calibre.lst"))
  calibre <- mutate(calibre, Basename = basename(filename),
                    ScanType = str_replace(Basename, r"{^T2_([^[:digit:]]+)_[[:digit:]].+\.nii.gz}", r"{\1}"))
  calibre <- mutate(calibre, PSFIm = file.path("calibre_psf_stepB", paste0("psf_", Basename)))
  #calibre <- mutate(calibre, PSFIm = file.path("calibre_psf", paste0("psf_", Basename)))
  #  calibre <- mutate(calibre, PSFIm = file.path("calibre_psf_uw", paste0("psf_", Basename)))
  
  # Load and extract psf profiles
  computePSF <- function(im, scantype) {
    img <- img_data(readNIfTI(im))
    extractPSF(img, scantype)
  }  
  computePSFIm <- function(im, scantype) {
    img <- img_data(readNIfTI(im))
    extractPSFIm(img, scantype)
  }  
  
  calibrePSF <- mutate(calibre, PSF = map2(PSFIm, ScanType, ~computePSF(.x, .y)))
  calibrePSFIm <- mutate(calibre, PSF = map2(PSFIm, ScanType, ~computePSFIm(.x, .y)))
  save(calibre, calibrePSF, calibrePSFIm, file=fl)
}

cdf <- function(x, mu, sigma){
  g <- dnorm(x, mean=mu, sd=sigma)
  return(g/max(g))
}
# 3d version
# Note that the manual for dmvnorm states that the sigma
# variable is a covariance matrix, while dnorm talks about standard deviations
# Sanity check on these two below:

library(minpack.lm)
library(mvtnorm)
cdf3 <- function(x, y, z, m1,m2,m3, s00,s01,s11,s02,s12,s22){
  xx <- as.matrix(cbind(x,y,z))
  mu <- c(m1,m2,m3)
  sigma <- matrix(c(s00,s01,s02,s01,s11,s12,s02,s12,s22), nrow=3)
  g <- dmvnorm(xx, mean=mu, sigma=sigma)
  return(g/max(g))
}

cdf3P <- function(x, y, z, m1,m2,m3, s00,s11,s22){
  # assuming PSF is axis aligned
  xx <- as.matrix(cbind(x,y,z))
  mu <- c(m1,m2,m3)
  sigma <- matrix(c(s00,0,0,0,s11,0,0,0,s22), nrow=3)
  g <- dmvnorm(xx, mean=mu, sigma=sigma)
  return(g/max(g))
}


cdf3PB <- function(x, y, z, m1,m2,m3, s00,s22){
  # assuming PSF is axis aligned and the in-plane sigmas
  # are the same
  xx <- as.matrix(cbind(x,y,z))
  mu <- c(m1,m2,m3)
  sigma <- matrix(c(s00,0,0,0,s00,0,0,0,s22), nrow=3)
  g <- dmvnorm(xx, mean=mu, sigma=sigma)
  return(g/max(g))
}

# This seems to failing for the "perfect" recommended data on some machines
# There's mention of this in the docs. Not sure why it works on n183.
# ignore for now - when it did work the results matched perfectly
#
# aa <- summarise(group_by(calibreProfilesFast, ScanType, direction, type), SD=sd(value), .groups="drop")
# aa <- mutate(aa, FWHM=SD*2*sqrt(2*log(2)))

doFit <- function(df) {
  fit <- nls(value ~ cdf(xcoord, mu, sigma),  
             start=list(mu=0, sigma=3), 
             lower=c(-1,0.01),
             upper=c(1, 10),
             algorithm="port",
             data=df)
  cc <- coef(fit)
  return(tibble(mu=cc[1], sigma=cc[2]))
}

doFit3d <- function(df) {
  startparams <- list(m1=0, m2=0, m3=0, s00=1, s01=0, s11=1, s02=0, s12=0, s22=1)
  fit <- nlsLM(brightness ~ cdf3(inplaneA, inplaneB, thruplane, m1, m2, m3, s00,s01,s11,s02,s12,s22),
               data=df, start=startparams)
  return(fit)
}

doFit3dP <- function(df) {
  # version parallel to axis
  # could also force the two inplane values to be the same
  startparams <- list(m1=0, m2=0, m3=0, s00=1, s11=1, s22=4.5)
  lower = c(-1,-1,-1, 0.5,0.5,0.5)
  upper = c(1, 1, 1, 3,3,10)
  fit <- nlsLM(brightness ~ cdf3P(inplaneA, inplaneB, thruplane, m1, m2, m3, s00,s11,s22),
               data=df, start=startparams, lower=lower, upper=upper)
  return(fit)
}
doFit3dPB <- function(df) {
  # version parallel to axis
  # could also force the two inplane values to be the same
  startparams <- list(m1=0, m2=0, m3=0, s00=1,  s22=4.5)
  lower = c(-1,-1,-1, 0.5,0.5)
  upper = c(1, 1, 1,3,10)
  fit <- nlsLM(brightness ~ cdf3PB(inplaneA, inplaneB, thruplane, m1, m2, m3, s00,s22),
               data=df, start=startparams, lower=lower, upper=upper)
  return(fit)
}

# each set of profiles is packed into 1 dataframe, unnest and then repack into smaller groups
calibrePSF.unnest <- unnest(select(calibrePSF, -ScanType), cols=c(PSF))
calibrePSF.repacked <- nest(group_by(calibrePSF.unnest, filename, direction, type, ScanType))

calibrePSF.repacked <- mutate(calibrePSF.repacked, gaussian_params = map(data, doFit))
calibreParams <- unnest(select(calibrePSF.repacked, -data), gaussian_params)
calibreParams <- separate_wider_delim(calibreParams, cols = ScanType, delim="_", too_few="align_start", names=c("orientation", "style"), cols_remove=FALSE)
calibreParams <- mutate(calibreParams, style=if_else(is.na(style), "Orig", style))
calibreParams <- mutate(calibreParams, FWHM =  sigma * 2*sqrt(2*log(2)))
#ggplot(calibreParams, aes(x=direction, y=sigma, colour=type))  + geom_point() + facet_wrap(~ScanType)

ca1 <- calibrePSF.unnest
ggplot(filter(ca1, type=="observed"), aes(x=xcoord, y=value, group = interaction(filename, direction, type, ScanType))) + 
 geom_line(aes(colour=ScanType)) + facet_grid(ScanType ~ direction)

ggplot(filter(calibrePSF.unnest, type=="observed"), aes(x=xcoord, y=value, group = interaction(filename, direction, type, ScanType))) + 
  geom_line(aes(colour=ScanType)) + facet_grid(ScanType ~ direction)

#ggplot(calibreParams, aes(x=direction, y=sigma, colour=type))  + geom_point() + facet_grid(style~orientation)
#ggplot(calibreParams, aes(x=direction, y=FWHM, colour=type))  + geom_point() + facet_grid(style~orientation)

#ggplot(calibreParams, aes(x=direction, y=FWHM, colour=type, group=interaction(type, direction)))  + geom_boxplot() + geom_point() + facet_grid(style~orientation)


# Lets assume that inplaneA and inplaneB are the same
calibreParams <- mutate(calibreParams, direction = if_else(direction == "through plane", direction, "inplane"))

#ggplot(calibreParams, aes(x=direction, y=FWHM, colour=type, group=interaction(type, direction)))  + geom_boxplot() + geom_point() + facet_grid(style~orientation)

calibreMeanParams <- summarise(ungroup(calibreParams), mFWHM=mean(FWHM), mSigma = mean(sigma), mVar = mean(sigma^2), .by=c(orientation, style, direction, type))

VoxelSpacing <- tibble(style = rep(c("Orig", "Fast"), 2), direction=rep(c("inplane", "through plane"), c(2,2)), spacing=c(1.5, 1.6, 5, 5))

calibreMeanParams <- left_join(ungroup(calibreMeanParams), VoxelSpacing)
calibreMeanParams <- mutate(calibreMeanParams, multiplier = mFWHM/spacing)

#View(filter(calibreMeanParams, type=="observed"))
#View(filter(calibreMeanParams, type=="observed", direction=="inplane"))
#View(filter(calibreMeanParams, type=="observed", direction=="through plane"))

# estimated parameters are quite different between scan types and orientations.
# Need support in the code for different PSFs

calibrePSFIm <- mutate(calibrePSFIm, params = map(PSF, doFit3dPB))
#calibrePSFIm <- mutate(calibrePSFIm, params = map(PSF, doFit3dP))
# These are surprisingly variable.
# Is registration failing?
# row 21 has a massive thru plane sigma
ggplot(filter(calibrePSFIm$PSF[[21]], thruplane==0), aes(x=inplaneA, y=inplaneB, fill=brightness)) + geom_tile()
ggplot(filter(calibrePSFIm$PSF[[21]], inplaneA==0), aes(x=inplaneB, y=thruplane, fill=brightness)) + geom_tile()

# registration wasn't perfect - I've tweaked it to improve things. It is better, but still not perfect.
# The distortion will limit the possible quality improvement
calibrePSFIm <- mutate(calibrePSFIm, sigmas = map(params, ~coef(.x)[4:6]))
calibrePSFIm <- mutate(calibrePSFIm, s00 = map_dbl(sigmas, "s00"), s22=map_dbl(sigmas, "s22"))


#calibrePSFIm <- mutate(calibrePSFIm, s00 = map_dbl(sigmas, "s00"), s11=map_dbl(sigmas, "s11"), s22=map_dbl(sigmas, "s22"))


summarise(calibrePSFIm, s00=mean(s00), s22=mean(s22), .by=ScanType)
summarise(calibrePSFIm, s00=mean(s00), s11=mean(s11), s22=mean(s22), .by=ScanType)

cc <- calibrePSFIm # stage 1

cc <- mutate(cc, Stage=1)
calibrePSFIm <- mutate(calibrePSFIm, Stage=2)
xx <- bind_rows(cc, calibrePSFIm)
cf <- filter(xx, ScanType=="COR_Fast")
ggplot(filter(cc$PSF[[1]], thruplane==0), aes(x=inplaneA, y=inplaneB, fill=brightness)) + geom_tile()
ggplot(filter(cc$PSF[[1]], inplaneA==0), aes(x=inplaneA, y=inplaneB, fill=brightness)) + geom_tile()

# Ignore this set - historical
# A tibble: 6 × 3
# ScanType   s00   s22
# <chr>    <dbl> <dbl>
#   1 AXI      0.713  3.60
#   2 AXI_Fast 0.852  3.02
#   3 COR      0.690  2.31
#   4 COR_Fast 1.07   2.54
#   5 SAG      0.976  3.25
#   6 SAG_Fast 1.17   3.60

# result after flipping the digital model in the inferior/superior direction
# and only including the copper solution in the model
# Pass 1
# 1 AXI      0.549  4.31
# 2 AXI_Fast 0.678  3.81
# 3 COR      0.5    3.37
# 4 COR_Fast 0.686  3.11
# 5 SAG      0.5    4.24
# 6 SAG_Fast 0.734  4.93

# Pass 2 - very similar, but coronal drops.
# I think this is because there are horizontal supports
# in the coronal plane that are visible, and seem to be
# biasing something. Use what we have
# ScanType   s00   s22
# <chr>    <dbl> <dbl>
# 1 AXI      0.540  4.25
# 2 AXI_Fast 0.700  3.43
# 3 COR      0.5    2.86
# 4 COR_Fast 0.694  2.62  # smaller again - why?
# 5 SAG      0.5    3.99
# 6 SAG_Fast 0.771  5.19

# A tibble: 6 × 3
# ScanType   s00   s22
# <chr>    <dbl> <dbl>
#   1 AXI      0.540  4.27
# 2 AXI_Fast 0.699  3.45
# 3 COR      0.625  1.64
# 4 COR_Fast 1.03   1.50
# 5 SAG      0.973  2.01
# 6 SAG_Fast 1.28   2.59


# This is the result from unwarping first - much larger - coronal
# isn't as different to the others
# ScanType   s00   s22
# <chr>    <dbl> <dbl>
#   1 AXI       1.96 10   
# 2 AXI_Fast  2.30 10   
# 3 COR       2.95  9.41
# 4 COR_Fast  3     7.95
# 5 SAG       2.12 10   
# 6 SAG_Fast  2.87 10   

# Strangely different in-plane values for cor and sag. Check their registration
# I've improved it as much as I can. I will force the optimizer to
# estimate a single in-plane sigma

# Sanity check - does dmvtnorm interpret sigma as covariances?
xx <- seq(from=-10, to=10, by=0.01)
yy <- dnorm(xx, sd=2)
yy <- yy/max(yy)
hh <- as.matrix(cbind(xx, 0))
# sigma of 2 corresponds to a cov of 4
zz <- dmvnorm(hh, sigma=4*diag(2))
zz <- zz/max(zz)
# OK - same, so it does take covariances as input.
range(yy-zz)
# perhaps the problem is driven by the ends - where the distortion is quite bad.
# Try generating results from the middle.
