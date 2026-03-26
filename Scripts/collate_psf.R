library(tidyverse)
library(oro.nifti)

# sag_fast - 5x1.6x1.6
# sag - 5x1.5x1.5
# cor_fast - 1.6x5x1.6
# cor - 1.5x5.1.5
# axi_fast - 1.6x1.6x5
# axi - 1.5x1.5x5
# focus on the data from calibre
fl <- "calibre_psf.Rda"
if (file.exists(fl)) {
  load(fl)
} else {
  calibre <- tibble(filename=readLines("calibre.lst"))
  calibre <- mutate(calibre, Basename = basename(filename),
                    ScanType = str_replace(Basename, r"{^T2_([^[:digit:]]+)_[[:digit:]].+\.nii.gz}", r"{\1}"))
  calibre <- mutate(calibre, PSFIm = file.path("calibre", paste0("psf_", Basename)))
  calibre <- nest(group_by(calibre, ScanType))
  # Load and average psf
  
  computePSF <- function(Ims) {
    ll <- length(Ims)
    total <- img_data(readNIfTI(Ims[1]))
    for (i in 2:ll) {
      total <- total + img_data(readNIfTI(Ims[i]))
    }
    total <- total/ll
    return(total)
  }
  calibrePSF <- mutate(calibre, meanPSF = map(data, ~computePSF(.x$PSFIm)))
  save(calibre, calibrePSF, file=fl)
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

calibrePSF <- mutate(ungroup(calibrePSF), Profiles = map2(meanPSF, ScanType, ~extractPSF(.x, .y)))

calibreProfiles <- bind_rows(calibrePSF$Profiles)

ggplot(calibreProfiles, aes(x=xcoord, y=value, colour=ScanType,group=interaction(direction,type, ScanType))) + geom_line() + facet_wrap(~direction)

ggplot(calibreProfiles, aes(x=xcoord, y=value, colour=ScanType,group=interaction(direction,type, ScanType))) + geom_line() + facet_grid(ScanType~direction)

ggplot(filter(calibreProfiles, type=="observed"), aes(x=xcoord, y=value, colour=ScanType,group=interaction(direction,type, ScanType))) + 
  geom_line() + 
  geom_line(data=filter(calibreProfiles, type=="recommended"), colour="black", alpha=0.5) +
  facet_grid(ScanType~direction)

calibreProfilesFast <- filter(calibreProfiles, str_detect(ScanType, "Fast"))
ggplot(filter(calibreProfilesFast, type=="observed"), aes(x=xcoord, y=value, colour=ScanType,group=interaction(direction,type, ScanType))) + 
  geom_line() + 
  geom_line(data=filter(calibreProfilesFast, type=="recommended"), colour="black", alpha=0.5) +
  facet_grid(ScanType~direction)

calibreProfilesOrig <- filter(calibreProfiles, !str_detect(ScanType, "Fast"))
ggplot(filter(calibreProfilesOrig, type=="observed"), aes(x=xcoord, y=value, colour=ScanType,group=interaction(direction,type, ScanType))) + 
  geom_line() + 
  geom_line(data=filter(calibreProfilesOrig, type=="recommended"), colour="black", alpha=0.5) +
  facet_grid(ScanType~direction)

# Looks like somewhat different results for original and fast sequences

# Original
# Through plane axial for original appears different to the others - about 1.5 times voxel size
# inplane - might be somewhat under done all round

# Fast
# Through plane looks OK
# in-planes needs to be broader - maybe about 1.5 instead of 1.2

# Estimating gaussian params for these
# 
cdf <- function(x, mu, sigma){
  g <- dnorm(x, mean=mu, sd=sigma)
  return(g/max(g))
}

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

calibreProfiles.nest <- nest(group_by(calibreProfiles, direction, ScanType, type))
calibreProfiles.nest <- mutate(calibreProfiles.nest, 
                               gaussian_params = map(data, doFit)
)


calibreParams <- unnest(select(calibreProfiles.nest, -data), gaussian_params)
calibreParams <- mutate(calibreParams, FWHM = sigma * 2*sqrt(2*log(2)))
calibreParams <- separate_wider_delim(calibreParams, cols = ScanType, delim="_", too_few="align_start", names=c("orientation", "style"), cols_remove=FALSE)
calibreParams <- mutate(calibreParams, style=if_else(is.na(style), "Orig", style))
ggplot(calibreParams, aes(x=direction, y=sigma, colour=type))  + geom_point() + facet_wrap(~ScanType)
ggplot(calibreParams, aes(x=direction, y=sigma, colour=type))  + geom_point() + facet_grid(style~orientation)

