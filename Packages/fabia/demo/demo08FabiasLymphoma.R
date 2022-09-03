#Data from Broad Institute "Cancer Program Data Sets"
#
#http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
#
#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"DLBCL_B"
#
#Rosenwald A, Wright G, Chan WC, Connors JM, Campo E, et al. (2002)
#  The use of molecular profiling to predict survival after chemotherapy
#  for diffuse large-B-cell lymphoma.
#N Engl J Med 346: 1937-1947.
#
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

data(DLBCL_B)

X <- as.matrix(XDLBCL)

resDLBCL2 <- fabias(X,5,0.6,300)

message("\n\nPlot of:
  1) data
  2) reconstructed data
  3) error (data - rec. data)
  4) absolute loadings
  5) absolute factors")
extractPlot(resDLBCL2,ti="FABIAS Lymphoma(Rosenwald)")

raDLBCL2 <- extractBic(resDLBCL2)

if ((raDLBCL2$bic[[1]][1]>1) && (raDLBCL2$bic[[1]][2])>1) {
message("\n\nPlot bicluster 1:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raDLBCL2,1)
}
if ((raDLBCL2$bic[[2]][1]>1) && (raDLBCL2$bic[[2]][2])>1) {
message("\n\nPlot bicluster 2:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raDLBCL2,2)
}

message("\n\nPlot of pairs of biclusters as biplots:
  1) rectangles are samples
     colors correspond to subclasses:
       oxidative phosphorylation, B-cell response, host response
     desired: to identify and separate the subclasses
  2) circles are genes
     red circles correspond to most indicative genes")
message("\n Plot1: Biclusters 1 and 2")
devAskNewPage(ask = TRUE)
plot(resDLBCL2,dim=c(1,2),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
message("\n Plot1: Biclusters 1 and 3")
plot(resDLBCL2,dim=c(1,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
message("\n Plot1: Biclusters 2 and 3")
plot(resDLBCL2,dim=c(2,3),label.tol=0.03,col.group=CDLBCL,lab.size=0.6)
devAskNewPage(ask = FALSE)
}

