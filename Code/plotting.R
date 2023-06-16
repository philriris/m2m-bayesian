library(rworldmap)
library(sf)
library(sp)
library(rgdal)
library(ggsn)
library(latex2exp)
library(patchwork)

load("vis_prep.RData")

dataset <- read.csv("locations.csv", header=TRUE)

base=getMap(resolution="low")
locations = unique(data.frame(SiteID = dataset$Date.n, 
                              Longitude = dataset$Longitude, 
                              Latitude = dataset$Latitude,
                              Where = dataset$River.basin))

rownames(locations) = as.character(locations$SiteID) 
locations= locations[,-1]
coordinates(locations) <-c("Longitude", "Latitude")
proj4string(locations) <- CRS("+proj=longlat +datum=WGS84")

base <- st_as_sf(base)
loc.sf <- st_as_sf(locations)
loc.bb <- st_as_sfc(st_bbox(loc.sf))

basin <- readOGR(dsn=".", layer = "amazon_basin") 
basin <- st_as_sf(basin)
rivers <- readOGR(dsn=".", layer = "amazon_rivers") 
rivers <- st_as_sf(rivers)

ylims <- c(st_bbox(loc.sf)[2]-2, st_bbox(loc.sf)[4]+2)

box <- st_bbox(loc.sf)
box[4] <- box[4]+2
box[2] <- box[2]-2
box[3] <- box[3]+1.5

main <- ggplot() + 
  geom_sf(data=base) +
  xlim(st_bbox(loc.sf)[c(1, 3)]) +
  ylim(ylims) +
  geom_sf(data=rivers, col="skyblue1") +
  geom_sf(data=loc.sf, col="grey30") + 
  geom_sf(data=base, col="black", fill=NA) +
  theme_bw() +
  theme(axis.title.x=element_blank(), 
        axis.title.y=element_blank()) +
  north(symbol = 10,
        x.min = box[1], x.max = box[3], y.min = box[2], y.max = box[4])

inset <- ggplot() + 
  geom_sf(data = base %>% dplyr::filter(GEO3 == "South America"), size=0.1, fill= "white") +
  geom_sf(data = basin, fill="palegreen2", color=NA, size=1.2) +
  geom_sf(data = base %>% dplyr::filter(GEO3 == "South America"), size=0.1, fill= NA) +
  geom_sf(data=st_as_sfc(box), fill=NA, col="black", size=1) +
  xlim(81,35) + theme_void() 

spdens <- ggplot(data=poly.spd$grid, aes(x = BPtoBCAD(calBP), y = PrDens)) +
  geom_area(fill="grey70") +
  geom_line(data = smooth.spd$grid, (aes(x = BPtoBCAD(calBP), y = PrDens))) +
  xlab("BC/AD") +
  ylab("Probability density") +
  theme_bw() +
  ggtitle("Summed probability of calibrated Polychrome Tradition dates")

sacomb <- 
  cowplot::ggdraw() +
  cowplot::draw_plot(main) +
  cowplot::draw_plot(inset, x = 0.75, y = 0.10, width = 0.3, height = 0.3)

polychrome_map <- sacomb / spdens +
  plot_layout(heights = c(6,2))

ggsave(plot = polychrome_map, file = "polychrome_map.pdf", width = 8.5, height = 8)

#SPD and KDE

pdf(file = "spd-ckde_polychrome.pdf", width = 10, height = 5) # save outputs

par(mfrow=c(1,2))

plot(poly.spd, calendar="BCAD")
plot(smooth.spd, calendar="BCAD", add=TRUE, type="simple")
plot(ckdepoly, calendar="BCAD", xlim=c(450,1850), type="multiline")

dev.off()

## Permutation testing

load("perm_res.RData")

pdf(file = "permtest_polychrome.pdf", width = 8, height = 4) # save outputs

par(mfrow=c(2,2))
par(mar=c(2,2.5,1.5,2))

plot(pt1, focalm="1", main= "Polychrome", calendar="BCAD")
plot(pt2, focalm="1", main= "Polychrome", calendar="BCAD")
plot(pt1, focalm="2", main= "Arauquinoid", calendar="BCAD")
plot(pt2, focalm="3", main= "Incised Rim", calendar="BCAD")

dev.off()

# Posterior check

load("mcmc_res.RData")

pdf(file = "posteriorcheck_polychrome.pdf", width = 7, height = 5) # save outputs
layout(matrix(c(
  1, 1, 1,
  2, 3, 4
), nrow=2, byrow=TRUE))

plot(pp.check, calendar="BCAD")

postHPDplot(chain_output[[1]]$chain1[,'r1'],xlab='',ylab='',
            show.hpd.val = FALSE,axes=F,
            hpd.col = "#F6A83A")
axis(side=1,cex.axis=0.5,padj=-1);mtext(TeX('$r1$'),side = 1,line=1.5,cex = 0.5)
abline(v=median(chain_output[[1]]$chain1[,"r1"]), lty=2, lwd=1.5, col="white")


postHPDplot(chain_output[[1]]$chain1[,'r2'],xlab='',ylab='',
            show.hpd.val = FALSE,axes=F,
            hpd.col = "#F6A83A")
axis(side=1,cex.axis=0.5,padj=-1);mtext(TeX('$r2$'),side = 1,line=1.5,cex = 0.5)
abline(v=median(chain_output[[1]]$chain1[,"r2"]), lty=2, lwd=1.5, col="white")


postHPDplot(chain_output[[1]]$chain1[,'chp'],xlab='',ylab='',
            show.hpd.val = FALSE,axes=F,
            hpd.col = "#BA9BC9")
axis(labels = paste(seq(500,700,50), "cal BP", sep=" "), at = seq(500,700,50),
     side=1,cex.axis=0.5,padj=-1);
labels <- paste("AD", BCADtoBP(seq(500,700,50)), sep=" ")
axis(labels = labels, at = seq(500,700,50),
     side=1,cex.axis=0.5,padj=0.5)
mtext(TeX('$Changepoint$'),side = 1,line=1.5,cex = 0.7,padj=1.1)
abline(v=median(chain_output[[1]]$chain1[,"chp"]), lty=2, lwd=1.5, col="white")

dev.off()

# MCMC diagnostics

pdf(file = "mcmc_diags.pdf", width = 8, height = 9) # save outputs

plot(chain_output[[1]], col=c("#F6A83A", "#BA9BC9", "palegreen2"))

dev.off()