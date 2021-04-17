library(phytools) ## should be phytools >= 0.7-76 https://github.com/liamrevell/phytools

## function to compute geometric mean
gmean<-function(x) exp(mean(log(x)))

## Figure 1
## a) Simple phylogeny of five taxa with the total edge length (above each branch) and 
## sigma^2 rates (below it). b) Calculation of the expected variance covariance matrix of tip 
## values for the trait x, T, given the branch length and rates of panel a). Each i,jth element 
## corresponds to the sum of the products above and below each edge along all the branches leading 
## from the root to the common ancestor of taxa i and j.
library(phytools)
set.seed(88)
tree<-pbtree(n=5,scale=10,tip.label=LETTERS[5:1])
sig2x<-exp(fastBM(tree,internal=TRUE))
tt<-tree
tt$edge.length<-tt$edge.length*apply(tt$edge,1,
    function(e,x) phytools:::ln.mean(x[e]),x=sig2x)
pdf(file="Figure1.pdf",width=10,height=6)
## layout(matrix(c(1,2),2,1),heights=c(0.6,0.4))
par(mfrow=c(1,2))
plotTree(tree,mar=c(3.1,1.1,4.1,1.1),fsize=1.25,
	ylim=c(0.5,5.2),lwd=7,color="black")
par(fg="transparent")
plotTree(tree,mar=c(3.1,1.1,4.1,1.1),fsize=1.25,
	ylim=c(0.5,5.2),lwd=5,color="grey",add=TRUE)
par(fg="black")
axis(1)
edgelabels(round(tree$edge.length,2),pos=3,
	frame="none",cex=0.7)
edgelabels(round(apply(tt$edge,1,
	function(e,x) phytools:::ln.mean(x[e]),x=sig2x),2),
	pos=1,frame="none",cex=0.7)
mtext("a)",line=1,adj=0)
plot.new()
par(mar=c(3.1,1.1,4.1,1.1))
plot.window(xlim=c(0,6),ylim=c(-4,6))
lines(c(0,6,6,0,0),c(0,0,6,6,0))
for(i in 1:5) lines(c(i,i),c(0,6))
for(i in 1:5) lines(c(0,6),c(i,i))
V<-vcv(tt)[LETTERS[1:5],LETTERS[1:5]]
for(i in 1:5) text(i+0.5,5.5,LETTERS[i],cex=1.1)
for(i in 1:5) text(0.5,5.5-i,LETTERS[i],cex=1.1)
for(i in 1:5) for(j in 1:5) text(0.5+i,5.5-j,
	round(V[i,j],2),cex=1.1)
mtext("b)",line=1,adj=0)
dev.off()
## end Figure 1

## Figure 2
## a) Simulated evolutionary rates, sigma^2, in which the logarithm of the rate of evolution 
## evolves by Brownian motion on the tree. b) The phylogeny of panel a) projected into a space 
## defined by time (on the horizontal axis) and a simulated trait vector, x, obtained using the 
## rates of panel a).
library(phytools)
set.seed(77)
tree<-pbtree(n=60)
sig2x<-exp(fastBM(tree,internal=TRUE))
tt<-tree
tt$edge.length<-tt$edge.length*apply(tt$edge,1,
    function(e,x) phytools:::ln.mean(x[e]),x=sig2x)
x<-fastBM(tt,internal=FALSE)
object<-list()
class(object)<-"multirateBM"
object$sig2<-sig2x
object$tree<-tree
object<-plot(object,plot=FALSE)
pdf(file="Figure2.pdf",width=12,height=8)
par(mfrow=c(1,2))
plot(object,ftype="i",fsize=0.7,mar=c(3.1,2.1,3.1,1.1),
	outline=TRUE,digits=3)
axis(1,at=round(seq(0,max(nodeHeights(tree)),length.out=5),2))
mtext("a)",at=-1)
par(mar=c(3.1,4.1,3.1,1.1))
phenogram(object$tree,x,colors=object$cols,
	spread.cost=c(1,0),fsize=0.7,ftype="i",
	las=1,xlab="")
mtext("b)",at=-0.5)
dev.off()
## end Figure 2

## Figure 3
## Estimated values of sigma^2 (y) compared to their known true values (x) for the data of Figure 
## 2 obtained using the penalized-likelihood approach with four different values of the penalty 
## term coefficient, lambda: a) lambda = 0.01; b) lambda = 0.1; c) lambda = 1; and d) lambda = 10. 
## The 1:1 line is indicated.
x<-x[tree$tip.label]
lambda<-c(0.01,0.1,1.0,10.0)
fits.x<-list()
for(i in 1:length(lambda)) 
	fits.x[[i]]<-multirateBM(tree,x,lambda=lambda[i])
pdf(file="Figure3.pdf",width=8,height=8)
par(mfrow=c(2,2))
for(i in 1:length(fits.x)){
	xylim<-range(c(sig2x,fits.x[[i]]$sig2))
	plot(sig2x,fits.x[[i]]$sig2,log="xy",xlim=xylim,,
		ylim=xylim,bty="n",las=1,pch=21,cex=1.2,
		bg="grey",xlab=expression(paste("true ",sigma^2)),
		ylab=expression(paste("estimated ",sigma^2)),
		cex.axis=0.8)
	lines(xylim,xylim)
	if(i==1) mtext(expression(paste("a) ",lambda,"= 0.01")),adj=0)
	else if(i==2) mtext(expression(paste("b) ",lambda,"= 0.1")),adj=0)
	else if(i==3) mtext(expression(paste("c) ",lambda,"= 1.0")),adj=0)
	else if(i==4) mtext(expression(paste("d) ",lambda,"= 10")),adj=0)
}
dev.off()
## end Figure 3

## Figure 4
## a) Simulated evolutionary rates, sigma^2, in which the logarithm of the rate of evolution 
## is uncorrelated between nodes and tips on the phylogeny. b) The tree projected into a space 
## defined by time (on the horizontal axis) and a simulated trait vector, x, obtained using the 
## rates of panel a).
library(phytools)
set.seed(77)
tree<-pbtree(n=60)
set.seed(66)
sig2y<-exp(rnorm(n=Ntip(tree)+tree$Nnode,sd=2.5))
tt<-tree
tt$edge.length<-tt$edge.length*apply(tt$edge,1,
    function(e,x) phytools::ln.mean(x[e]),x=sig2y)
y<-fastBM(tt)
object<-list()
class(object)<-"multirateBM"
object$sig2<-sig2y
object$tree<-tree
object<-plot(object,plot=FALSE,digits=3)
pdf(file="Figure4.pdf",width=12,height=8)
par(mfrow=c(1,2))
plot(object,ftype="i",fsize=0.7,mar=c(3.1,2.1,3.1,1.1),
	outline=TRUE,digits=3)
axis(1,at=round(seq(0,max(nodeHeights(tree)),length.out=5),2))
mtext("a)",at=-1)
par(mar=c(3.1,4.1,3.1,1.1))
phenogram(object$tree,y,colors=object$cols,
	spread.cost=c(1,0),fsize=0.7,ftype="i",
	las=1,xlab="")
mtext("b)",at=-0.5)
dev.off()
## end Figure 4

## Figure 5
## Estimated values of sigma^2 (y) compared to their known true values (x) for the data of 
## Figure 4 obtained using the penalized-likelihood approach using four different values of the 
## penalty term coefficient, lambda: a) lambda = 0.01; b) lambda = 0.1; c) lambda = 1; and d) 
## lambda = 10. The 1:1 line is indicated.
lambda<-c(0.01,0.1,1.0,10.0)
fits.y<-list()
for(i in 1:length(lambda)) 
	fits.y[[i]]<-multirateBM(tree,y,lambda=lambda[i])
pdf(file="Figure5.pdf",width=8,height=8)
par(mfrow=c(2,2))
for(i in 1:length(fits.y)){
	xylim<-range(c(sig2y,fits.y[[i]]$sig2))
	plot(sig2y,fits.y[[i]]$sig2,log="xy",xlim=xylim,,
		ylim=xylim,bty="n",las=1,pch=21,cex=1.2,
		bg="grey",xlab=expression(paste("true ",sigma^2)),
		ylab=expression(paste("estimated ",sigma^2)),
		cex.axis=0.8)
	lines(xylim,xylim)
	if(i==1) mtext(expression(paste("a) ",lambda,"= 0.01")),adj=0)
	else if(i==2) mtext(expression(paste("b) ",lambda,"= 0.1")),adj=0)
	else if(i==3) mtext(expression(paste("c) ",lambda,"= 1.0")),adj=0)
	else if(i==4) mtext(expression(paste("d) ",lambda,"= 10")),adj=0)
}
dev.off()
## end Figure 5

## Figure 6
## a) Simulated evolutionary rates, sigma^2, in which rate of evolution shifts discretely under 
## a continuous-time Markov process. b) The tree of panel a) projected into a space defined by time 
## (on the horizontal axis) and a simulated trait vector, x, obtained using the rates of panel a).
library(phytools)
set.seed(77)
tree<-pbtree(n=60)
set.seed(55)
Q<-matrix(c(-0.2,0.2,0,
	0.2,-0.4,0.2,
	0,0.2,-0.2),3,3,
	byrow=TRUE,
	dimnames=list(1:3,1:3))
map<-sim.history(tree,Q,anc="2",quiet=TRUE)
while(any(sapply(map$maps,length)>2)) 
	map<-sim.history(tree,Q,anc="2",quiet=TRUE)
for(i in 1:length(map$maps)){
	if(length(map$maps[[i]])>1){
		map$maps[[i]]<-setNames(rep(map$edge.length[i]/2,2),
			names(map$maps[[i]]))
	}
}
map<-read.simmap(text=capture.output(write.simmap(map,file="")))
z<-sim.rates(map,sig2=setNames(c(0.1,1,10),1:3))
cols<-setNames(colorRampPalette(c("blue","red"))(3),1:3)
plot(map)
lastPP <- get("last_plot.phylo", envir = .PlotPhyloEnv)
xlim<-c(-0.3, 1.05) * diff(lastPP$x.lim)
pdf(file="Figure6.pdf",width=12,height=8)
par(mfrow=c(1,2))
plot(map,colors=cols,ftype="i",fsize=0.7,
	mar=c(3.1,2.1,3.1,1.1),
	outline=TRUE,xlim=xlim)
axis(1,at=round(seq(0,max(nodeHeights(tree)),length.out=5),2))
mtext("a)",at=-1)
legend("topleft",legend=c("","",
	expression(paste(sigma^2," = 0.1")),
	expression(paste(sigma^2," = 1.0")),
	expression(paste(sigma^2," = 10"))),
	pch=22,col=c("transparent","transparent","black","black","black"),
	pt.bg=c("transparent","transparent",cols),pt.cex=2,bty="n")
par(mar=c(3.1,4.1,3.1,1.1))
phenogram(map,z,colors=cols,
	spread.cost=c(1,0),fsize=0.7,ftype="i",
	xlab="")
dev.off()
## end Figure 6

## Figure 7
## Estimated values of sigma^2 (y) compared to their known true values (x) for the data of 
## Figure 6 obtained using the penalized-likelihood approach using four different values of the 
## penalty term coefficient, lambda: a) lambda = 0.01; b) lambda = 0.1; c) lambda = 1; and d) 
## lambda = 10. The 1:1 line is indicated.
sig2z<-c(getStates(map,"tips"),getStates(map,"nodes"))
sig2z<-setNames(setNames(c(0.1,1,10),1:3)[sig2z],names(sig2z))
lambda<-c(0.01,0.1,1.0,10.0)
fits.z<-list()
for(i in 1:length(lambda)) 
	fits.z[[i]]<-multirateBM(tree,z,lambda=lambda[i])
pdf(file="Figure7.pdf",width=8,height=8)
par(mfrow=c(2,2))
for(i in 1:length(fits.z)){
	xylim<-range(c(sig2z,fits.z[[i]]$sig2))
	plot(sig2z,fits.z[[i]]$sig2,log="xy",xlim=xylim,,
		ylim=xylim,bty="n",las=1,pch=16,cex=1.2,
		col=make.transparent("grey",0.5),
		xlab=expression(paste("true ",sigma^2)),
		ylab=expression(paste("estimated ",sigma^2)),
		cex.axis=0.8)
	levs<-c(0.1,1,10)
	lines(levs,sapply(levs,function(x,y,z) mean(y[z==x]),
		y=fits.z[[i]]$sig2,z=sig2z),type="b",pch=16,
		cex=1.5)
	lines(xylim,xylim)
	if(i==1) mtext(expression(paste("a) ",lambda,"= 0.01")),adj=0)
	else if(i==2) mtext(expression(paste("b) ",lambda,"= 0.1")),adj=0)
	else if(i==3) mtext(expression(paste("c) ",lambda,"= 1.0")),adj=0)
	else if(i==4) mtext(expression(paste("d) ",lambda,"= 10")),adj=0)
}
dev.off()
## end Figure 7

## Figure 8
## Estimated values of sigma^2 (y) compared to their known true values (x) for data simulated 
## with a constant rate of evolution, sigma^2 = 1. Vertical and horizontal lines show the true 
## rate of evolution from the simulation, and the Maximum Likelihood estimate of the rate obtained 
## in a single-rate analysis.
library(phytools)
set.seed(77)
tree<-pbtree(n=60)
set.seed(44)
zz<-fastBM(tree)
lambda<-c(0.01,0.1,1.0,10.0)
fits.zz<-list()
for(i in 1:length(lambda)) 
	fits.zz[[i]]<-multirateBM(tree,zz,lambda=lambda[i])
pdf(file="Figure8.pdf",width=8,height=8)
par(mfrow=c(2,2))
xylim<-range(c(1,sapply(fits.zz, function(x) x$sig2)))
for(i in 1:length(fits.zz)){
	plot(NA,log="xy",xlim=xylim,,
		ylim=xylim,bty="n",las=1,
		xlab=expression(paste("true ",sigma^2)),
		ylab=expression(paste("estimated ",sigma^2)),
		cex.axis=0.8)
	abline(v=1,lty="dotted")
	abline(h=geiger::fitContinuous(tree,zz)$opt$sigsq,
		lty="dotted")
	points(rep(1,length(fits.zz[[i]]$sig2)),
		fits.zz[[i]]$sig2,pch=16,cex=1.2,
		col=make.transparent("grey",0.5))
	lower<-gmean(fits.zz[[i]]$sig2)-sqrt(var(fits.zz[[i]]$sig2))
	upper<-gmean(fits.zz[[i]]$sig2)+sqrt(var(fits.zz[[i]]$sig2))
	points(1,mean(fits.zz[[i]]$sig2),type="b",pch=16,
		cex=1.5)
	if(i==1) mtext(expression(paste("a) ",lambda,"= 0.01")),adj=0)
	else if(i==2) mtext(expression(paste("b) ",lambda,"= 0.1")),adj=0)
	else if(i==3) mtext(expression(paste("c) ",lambda,"= 1.0")),adj=0)
	else if(i==4) mtext(expression(paste("d) ",lambda,"= 10")),adj=0)
}
dev.off()
## end Figure 8

## Figure 9
## Fitted variable-rate Brownian motion model for log(body mass) evolution in 49 species of 
## mammals. Different values of lambda correspond to different penalty coefficients for rate 
## variation among edges of the tree. Note that each panel has a different scale.
data(mammal.tree)
data(mammal.data)
lnSize<-setNames(log(mammal.data$bodyMass),rownames(mammal.data))
lambda<-c(0.01,0.1,1.0,10.0)
fits.mammals<-list()
for(i in 1:length(lambda))
	fits.mammals[[i]]<-multirateBM(mammal.tree,lnSize,lambda=lambda[i])
pdf(file="Figure9.pdf",width=10,height=10)
par(mfrow=c(2,2))
for(i in 1:length(fits.mammals)){
	object<-plot(fits.mammals[[i]],plot=FALSE)
	object<-setMap(object,heat.colors(n=1000)[1000:1])
	plot(object,fsize=0.7,outline=TRUE,mar=c(1.1,1.1,3.1,1.1),ftype="i")
	if(i==1) mtext(expression(paste("a) ",lambda,"= 0.01")),adj=0)
	else if(i==2) mtext(expression(paste("b) ",lambda,"= 0.1")),adj=0)
	else if(i==3) mtext(expression(paste("c) ",lambda,"= 1.0")),adj=0)
	else if(i==4) mtext(expression(paste("d) ",lambda,"= 10")),adj=0)
}
dev.off()
## end Figure 9

