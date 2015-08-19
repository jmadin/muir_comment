# Load original data sets downloaded from Dryad
spp <- read.csv("data/Valid species numbers.csv", as.is=TRUE)
lat <- read.csv("data/lat_by_depth_by_sp.csv", as.is=TRUE)

# Merge data sets
dat <- merge(spp, lat)

# Function for calculating the number of unique objects in a vector
length_unique <- function(x) {
  length(unique(x))
}

# The data contain two extra species, because only 104 in Muir et al. analysis
aggregate(dat$Latitude, list(dat$Species), length_unique)

# I. togianensis and A. russelli only at three latitudes and need to be removed

dat <- dat[dat$Species != "I. togianensis",]
dat <- dat[dat$Species != "A. russelli",]
length(unique(dat$Species))

# Center latitude at 35 to avoid dealing with negatives
dat$lat <- dat$Latitude + 35

# Create new latitude squared variable
dat$lat_sq <- dat$lat^2

# Load the package used by Muir et al. for random effect in a quantile regression model
library(lqmm)

# The basic Muir et al. model
rqm1 <- lqmm(Depth ~ lat + lat_sq, random = ~ 1, group=Species, tau=0.975, data=dat, control = lqmmControl(method = "df", UP_max_iter = 200))

# Quick check to see if their model best fixed-effects
rqm2 <- lqmm(Depth ~ lat, random = ~ 1, group=Species, tau=0.975, data=dat, control = lqmmControl(method = "df", UP_max_iter = 200))
rqm3 <- lqmm(Depth ~ 1, random = ~ 1, group=Species, tau=0.975, data=dat, control = lqmmControl(method = "df", UP_max_iter = 200))

AIC(rqm1) # 110324.7, best
AIC(rqm2) # 110816.9
AIC(rqm3) # 111020.4

# Check Muir et al. assumption about latitude-species random effects (see our comment for justification)

rqm4 <- lqmm(Depth ~ lat + lat_sq, random = ~ 1 + lat, group=Species, tau=0.975, data=dat, control = lqmmControl(method = "df", UP_max_iter = 200))
rqm5 <- lqmm(Depth ~ lat + lat_sq, random = ~ 1 + lat + lat_sq, group=Species, tau=0.975, data=dat, control = lqmmControl(method = "df", UP_max_iter = 200))

AIC(rqm4) # 109725
AIC(rqm5) # 109703.3, AIC is 621.4 better than rqm1

# Log-likelihood test

logLik(rqm5)
logLik(rqm1)

# Chi-sq
2*c(logLik(rqm4)-logLik(rqm1))

# p-value
pchisq(2*c(logLik(rqm5)-logLik(rqm1)), 7-5, lower.tail=FALSE)

# If species-latitude effects important, then look at them as fixed effects

library(quantreg)

# Some species only a few latitudes, so only species with n > 100
temp <- aggregate(dat$Depth, list(dat$Species), length)
spp_rm <- temp[temp$x < 100,]$Group.1
dat2<- dat
for (i in spp_rm) {
  dat2 <- dat2[dat2$Species != i,]
}

rq1 <- rq(Depth ~ lat:Species + lat_sq:Species + Species, tau=0.975, data=dat2)

AIC(rq1) # 93133.17

rq1_coefs <- summary(rq1)$coef # Gives warning: In summary.rq(rq1) : 1358 non-positive fis
1358/dim(dat2)[1] # ~11% non-positive fis relative to replicate is not considered large, see FAQ for quantreg

# Nonetheless, see if estimates are similar to running quantreg on each species separately

rq2_coefs <- c()
for (s in sort(unique(dat2$Species))) {
  rq_temp <- rq(Depth ~ lat + lat_sq, tau=0.975, data=dat2[dat2$Species == s,])
  rq2_coefs <- rbind(rq2_coefs, rq_temp$coef)
}

plot(rq2_coefs[2:50,1], as.vector(rq1_coefs[2:50,1]))
plot(rq2_coefs[,2], as.vector(rq1_coefs[51:100,1]))
plot(rq2_coefs[,3], as.vector(rq1_coefs[101:150,1]))

# Estimates scale the same regardless of if species analysed separately or in same rq analysis

# Tally concave up and down

sum(rq1_coefs[101:150,1] > 0 & rq1_coefs[101:150,4] < 0.05) # 4 significant quadratic and opposite of Muir
sum(rq1_coefs[101:150,1] < 0 & rq1_coefs[101:150,4] < 0.05) # 17 significnat quadratic and follow Muir
sum(rq1_coefs[101:150,4] > 0.05)                            # 29 no significant quadratic pattern

sum(rq1_coefs[101:150,1] < 0) # 38 Number following Muir et al. pattern
sum(rq1_coefs[101:150,1] > 0) # 12 (24%) number opposite pattern

# Make a figure of species effects centered at equator
pdf("figs/muir_species_rq.pdf", height=5, width=8)

  plot(dat2$lat-35, -dat2$Depth, ylim=c(-30, 30), xlim=c(-35, 35), col="grey", xlab="Latitude", ylab="Maximum depth - maximum depth at Equator (m)", axes=FALSE, type="n")
  axis(2, las=2)
  axis(1, at=seq(-30, 30, 30))

  # Add species lines
  for (i in unique(sort(dat2$Species))) {
    ss = seq(round(min(dat2$lat[dat2$Species==i])), round(max(dat2$lat[dat2$Species==i])), 0.5)
    pp <- predict(rq1, list(lat=ss, lat_sq=ss^2, Species=rep(i, length(ss))))
    # temp <- dat2[dat2$Species == i,]
    # rq_temp <- rq(Depth ~ lat + lat_sq, tau=0.975, data=temp)
    # pp <- predict(rq_temp, list(lat=ss, lat_sq=ss^2))
    ppp <- pp - pp[which(ss==35)]
    lines(ss-35, -ppp, col=rgb(0, 0, 0, 0.3))
  }

  # Add Muir et al. as thicker line
  ord <- order(dat$lat)
  pp <- predict(rqm1, 0)
  cor_lat <- dat$lat[ord]-35
  cor_dep <- -pp[ord]
  lines(cor_lat, cor_dep-cor_dep[which(cor_lat==0)], lwd=3, col="black")
  box()

dev.off()
