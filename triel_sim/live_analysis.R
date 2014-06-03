

particle0=read.csv("partInfo0.csv")
particle1=read.csv("partInfo1.csv")
particle2=read.csv("partInfo2.csv")

zoomX <- 0.8
zoomY <- 0.5

ct <- particle0[[13]][length(particle0[[13]])]
c0posx <- particle0[[1]][length(particle0[[1]])]
c1posx <- particle1[[1]][length(particle0[[1]])]
c2posx <- particle2[[1]][length(particle0[[1]])]

c0vx <- particle0[[4]][length(particle0[[4]])]
c1vx <- particle1[[4]][length(particle0[[4]])]
c2vx <- particle2[[4]][length(particle0[[4]])]



#interpolatedp2x <- seq.int(2, 6.9999, 0.0001)
par(mfrow=c(2,2))
plot(data.frame(particle0[13], particle0[1]), xlim=c(ct*zoomX, ct), ylim=c(c0posx*zoomY, c0posx+(c0posx*zoomY)))
points(data.frame(particle1[13], particle1[1]), col=2)
points(data.frame(particle2[13], particle2[1]), col=3)
legend("topleft",c("p0","p1","p2"), col=c(1,2,3), pch=c(1,1,1))


plot(data.frame(particle0[13], particle0[4]), xlim=c(ct*zoomX, ct), ylim=c(c0vx*zoomY, c0vx + (c0vx*zoomY)))
points(data.frame(particle1[13], particle1[4]), col=2)
points(data.frame(particle2[13], particle2[4]), col=3)
points(data.frame(particle2[13], particle2[10]/0.0010229858259037613), col=4)
legend("topleft",c("p0","p1","p2","p2lb"), col=c(1,2,3,4), pch=c(1,1,1,3))


l0 <- sqrt(((particle1[[1]] - particle0[[1]])^2) + (particle1[[2]] - particle0[[2]])^2)
lp0 <- sqrt(((particle2[[1]] - particle0[[1]])^2) + (particle2[[2]] - particle0[[2]])^2)
## #interpolatedlp0 <- sqrt(((interpolatedp2x - particle0[[posx]])^2) + (rep(5,50000) - particle0[[posy"]])^2) 

## ## l0 <- l0[!is.na(l0)]
## ## lp0 <- lp0[!is.na(lp0)]


plot(particle0[[13]], l0, xlim=c(ct*zoomX, ct))
points(particle0[[13]], lp0, col=2)
points(data.frame(particle0[13], (0.5 * l0 * lp0)), col=3)
legend("topleft",c("l(0)","l'(0)","A(0)"), col=c(1,2,3), pch=c(1,1,1))
