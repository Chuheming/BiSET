#
#
# Author: SEPP HOCHREITER
###############################################################################

n = 1000
l= 100
p = 10

dat <- makeFabiaDataBlocks(n = n,l= l,p = p,f1 = 5,f2 = 5,
  of1 = 5,of2 = 10,sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,
  sd_z = 1.0,sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

X <- dat[[1]]
Y <- dat[[2]]
ZC <- dat[[3]]
LC <- dat[[4]]

gclab <- rep.int(0,l)
gllab <- rep.int(0,n)
clab <- as.character(1:l)
llab <- as.character(1:n)
for (i in 1:p){
 for (j in ZC[i]){
     clab[j] <- paste(as.character(i),"_",clab[j],sep="")
 }
 for (j in LC[i]){
     llab[j] <- paste(as.character(i),"_",llab[j],sep="")
 }
 gclab[unlist(ZC[i])] <- gclab[unlist(ZC[i])] + p^i
 gllab[unlist(LC[i])] <- gllab[unlist(LC[i])] + p^i
}


groups <- gclab

#### FABIAS

resToy2 <- fabias(X,13,0.6,400)

message("\n\nPlot of:
  1) noise free data
  2) data
  3) reconstructed data
  4) error (data - rec. data)
  5) absolute loadings
  6) absolute factors")
extractPlot(resToy2,ti="FABIAS",Y=Y)


raToy2 <- extractBic(resToy2)

if ((raToy2$bic[[1]][1]>1) && (raToy2$bic[[1]][2])>1) {
message("\n\nPlot bicluster 1:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raToy2,1)
}
if ((raToy2$bic[[2]][1]>1) && (raToy2$bic[[2]][2])>1) {
message("\n\nPlot bicluster 2:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raToy2,2)
}

colnames(X(resToy2)) <- clab

rownames(X(resToy2)) <- llab


message("\n\nPlot of pairs of biclusters as biplots:
  1) rectangles are samples
     colors correspond to different biclusters membership
     labeled by: bicluster1_bicluster2_..._sampleID
  2) circles are genes
     indicative are large and red
     labeled by: bicluster1_bicluster2_..._geneID")
message("\n Plot1: Biclusters 1 and 2")
devAskNewPage(ask = TRUE)
plot(resToy2,dim=c(1,2),label.tol=10,col.group = groups,lab.size=0.6)
message("\n Plot1: Biclusters 1 and 3")
plot(resToy2,dim=c(1,3),label.tol=10,col.group = groups,lab.size=0.6)
message("\n Plot1: Biclusters 2 and 3")
plot(resToy2,dim=c(2,3),label.tol=10,col.group = groups,lab.size=0.6)

