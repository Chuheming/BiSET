#Data from Broad Institute "Cancer Program Data Sets"
#
#http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
#
#Hoshida Y, Brunet J-P, Tamayo P, Golub TR, Mesirov JP (2007)
#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"Breast_A" without Array S54
#van't Veer LJ, Dai H, van de Vijver MJ, He YD, Hart AA, et al. (2002)
#Gene expression profiling predicts clinical outcome of breast cancer.
#Nature 415:530-536.
#
# Author: SEPP HOCHREITER
###############################################################################

avail <- require(fabiaData)

if (!avail) {
    message("")
    message("")
    message("#####################################################")
    message("Package 'fabiaData' is not available: please install.")
    message("#####################################################")
} else {

data(Breast_A)

X <- as.matrix(XBreast)


resBreast1 <- fabia(X,5,0.1,400)

message("\n\nPlot of:
  1) data
  2) reconstructed data
  3) error (data - rec. data)
  4) absolute loadings
  5) absolute factors")
extractPlot(resBreast1,ti="FABIA Breast cancer(Veer)")


raBreast1 <- extractBic(resBreast1)

if ((raBreast1$bic[[1]][1]>1) && (raBreast1$bic[[1]][2])>1) {
message("\n\nPlot bicluster 1:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raBreast1,1)
}
if ((raBreast1$bic[[2]][1]>1) && (raBreast1$bic[[2]][2])>1) {
message("\n\nPlot bicluster 2:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raBreast1,2)
}


message("\n\nPlot of pairs of biclusters as biplots:
  1) rectangles are samples
     color corresponds to previously identified subgroups
     desired: separation of the subgroups
  2) circles are genes
     red circles correspond to most indicative genes")
message("\n Plot1: Biclusters 1 and 2")
devAskNewPage(ask = TRUE)
plot(resBreast1,dim=c(1,2),label.tol=0.03,col.group=CBreast,lab.size=0.6)
message("\n Plot1: Biclusters 1 and 3")
plot(resBreast1,dim=c(1,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
message("\n Plot1: Biclusters 2 and 3")
plot(resBreast1,dim=c(2,3),label.tol=0.03,col.group=CBreast,lab.size=0.6)
devAskNewPage(ask = FALSE)
}
