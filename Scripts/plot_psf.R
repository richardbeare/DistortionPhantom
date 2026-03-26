library(tidyverse)
library(oro.nifti)

psf_h <- oro.nifti::readNIfTI("mean_ax.nii.gz")
psf <- img_data(psf_h)

d <- dim(psf)

md <- d/2 + 1
bb <- 35
bbpo <- bb + 1
cc <- 2*bb + 1
begin <- md - bb
theend <- md + bb

middle <- psf[begin[1]:theend[1], begin[2]:theend[2], begin[3]:theend[3]]

# Says : in plane PSF is a sinc function with FWHM = 1.2xin plane res.
# Through-plane is gaussian with FWHM= slice thickness

# Note: FWHM = 2 * sqrt(2 * ln(2)) * sigma
# These are the same as used in NiftyMIC
sig_inplane <- (1.2 * 1.5)/(2*sqrt(2*log(2)))
sig_thruplane <- 5/(2*sqrt(2*log(2)))

gaussian <- function(X, sigma, mu=0) {
  1.0/(sigma*sqrt(2*pi)) * exp(-0.5 * (X - mu)^2/sigma^2)
}

dftheory <- tibble( xcoord = 0.5*((1:cc) - bb))
dftheoryA <- mutate(dftheory, value = gaussian(xcoord, sig_inplane), value=value/max(value), direction="inplaneA", type="recommended")
dftheoryB <- mutate(dftheory, value = gaussian(xcoord, sig_inplane), value=value/max(value), direction="inplaneB", type="recommended")
dftheoryC <- mutate(dftheory, value = gaussian(xcoord, sig_thruplane), value=value/max(value), direction="through plane", type="recommended")

dfA <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,bbpo,], direction="through plane", type="observed")
dfA <- mutate(dfA, value = value/max(value))
dfB <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[bbpo,,bbpo], direction="inplaneA", type="observed")
dfB <- mutate(dfB, value = value/max(value))

dfC <- tibble(xcoord = 0.5*((1:cc) - bb), value = middle[,bbpo,bbpo], direction="inplaneB", type="observed")
dfC <- mutate(dfC, value = value/max(value))

df <- bind_rows(dfA, dfB, dfC, dftheoryA, dftheoryB, dftheoryC)

ggplot(df, aes(x=xcoord, y=value, colour=type, group=interaction(direction, type))) + geom_point() + geom_line() +
  facet_wrap(~direction) +
  ylab("weight") + xlab("spread (mm)") + labs(colour=element_blank()) + ggtitle("Observed and recommended PSF")
ggsave("psf_plot.png")
