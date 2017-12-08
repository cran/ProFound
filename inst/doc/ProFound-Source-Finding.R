## ---- eval=FALSE---------------------------------------------------------
#  library(devtools)
#  install_github('asgr/ProFound')

## ------------------------------------------------------------------------
library(ProFound)

## ------------------------------------------------------------------------
image=readFITS(system.file("extdata", 'VIKING/mystery_VIKING_Z.fits',package="ProFound"))

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimageWCS(image)

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
out_segim=segim=profoundMakeSegim(image$imDat, magzero=30, pixscale=0.339, header=image$hdr, plot=TRUE)
out_profound=profoundProFound(image, magzero=30, verbose=TRUE, plot=TRUE)

## ------------------------------------------------------------------------
out_profound$segstats[1:10,]

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magplot(out_segim$segstats[1:40,c("R50", "SB_N90")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), log='x', xlim=c(0.4,5), ylim=c(22,25), grid=TRUE, xlab='R50 / asec', ylab='mag / asec^2')
points(out_profound$segstats[1:40,c("R50", "SB_N90")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), pch=16)
arrows(out_segim$segstats$R50[1:40], out_segim$segstats$SB_N90[1:40], out_profound$segstats$R50[1:40], out_profound$segstats$SB_N90[1:40], col='lightgrey', length=0)
rect(0.9, 23.5, 1.3, 24.3)
legend('bottomleft', legend=c('profoundProFound', 'profoundMakeSegim'), pch=c(16,1))
magbar('topright', title='Axrat', titleshift=1)

magplot(out_segim$segstats[1:40,c("R50", "con")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), log='x', xlim=c(0.4,5), ylim=c(0,1), grid=TRUE, xlab='R50 / asec', ylab='Concentration')
points(out_profound$segstats[1:40,c("R50", "con")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), pch=16)
arrows(out_segim$segstats$R50[1:40], out_segim$segstats$con[1:40], out_profound$segstats$R50[1:40], out_profound$segstats$con[1:40], col='lightgrey', length=0)
rect(0.9, 0.4, 1.3, 0.6)
legend('bottomleft', legend=c('profoundProFound', 'profoundMakeSegim'), pch=c(16,1))
magbar('topright', title='Axrat', titleshift=1)

magplot(out_segim$segstats[1:40,c("R50", "mag")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), log='x', xlim=c(0.4,5), ylim=c(17,24), grid=TRUE, xlab='R50 / asec', ylab='Mag')
points(out_profound$segstats[1:40,c("R50", "mag")], col=hsv(magmap(out_segim$segstats$axrat, flip=TRUE)$map), pch=16)
arrows(out_segim$segstats$R50[1:40], out_segim$segstats$mag[1:40], out_profound$segstats$R50[1:40], out_profound$segstats$mag[1:40], col='lightgrey', length=0)
rect(0.9, 20, 1.3, 22)
legend('bottomleft', legend=c('profoundProFound', 'profoundMakeSegim'), pch=c(16,1))
magbar('topright', title='Axrat', titleshift=1)

## ---- fig.width=6, fig.height=5------------------------------------------
maghist(out_profound$sky, xlab='Sky')
maghist(out_profound$skyRMS, xlab='Sky RMS')

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimageWCS(image)
magimageWCS(out_profound$sky, image$hdr)
magimageWCS(out_profound$skyRMS, image$hdr)

## ---- fig.width=5, fig.height=4------------------------------------------
maghist(out_profound$segstats[,'N100']/out_profound$segstats[,'flux'], xlab=' Worst case fraction error')

## ------------------------------------------------------------------------
out_profound=profoundProFound(image, magzero=30, verbose=TRUE, boundstats=TRUE)

## ------------------------------------------------------------------------
out_profound$segstats[1:10,c("Nedge", "Nsky", "Nobject", "Nborder", "edge_frac", "edge_excess", "flag_border")]

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_frac==1,"segID"]), col=c(0,rainbow(100)))
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_frac!=1,"segID"]), col=c(0,rainbow(100)))

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$Nborder==0,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Isolated Segment')
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$Nborder>0,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Non_Isolated Segment')

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_excess<1,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Fairly Elliptical Segment')
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_excess>1,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Non-Elliptical Segment')

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$flag_border==1,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Segment Borders Bottom')
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$flag_border==2,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Segment Borders Left')
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$flag_border==4,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Segment Borders Top')

## ---- fig.width=6, fig.height=6, dpi=40----------------------------------
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_frac>0.8,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Reliable Photometry Segments')
magimage(out_profound$segim*(out_profound$segim %in% out_profound$segstats[out_profound$segstats$edge_frac<0.8,"segID"]), col=c(0,rainbow(100)))
legend('topleft',legend='Unreliable Photometry Segments')

