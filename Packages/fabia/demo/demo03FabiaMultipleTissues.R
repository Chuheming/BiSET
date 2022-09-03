#Data from Broad Institute "Cancer Program Data Sets"
#
#http://www.broadinstitute.org/cgi-bin/cancer/datasets.cgi
#
#Subclass Mapping: Identifying Common Subtypes in Independent Disease Data Sets
#PLoS ONE 2(11): e1195. doi:10.1371/journal.pone.0001195
#
#"Multi_A" without Array BR_U16
#
#Su AI, Cooke MP, Ching KA, Hakak Y, Walker JR, et al. (2002)
#Large-scale analysis of the human and mouse transcriptomes.
#Proc Natl Acad Sci U S A 99:4465-4470.
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


data(Multi_A)

X <- as.matrix(XMulti)

resMulti1 <- fabia(X,5,0.06,300,norm=2)

message("\n\nPlot of:
  1) data
  2) reconstructed data
  3) error (data - rec. data)
  4) absolute loadings
  5) absolute factors")
extractPlot(resMulti1,ti="FABIA Multiple tissues(Su)")

raMulti1 <- extractBic(resMulti1)

if ((raMulti1$bic[[1]][1]>1) && (raMulti1$bic[[1]][2])>1) {
message("\n\nPlot bicluster 1:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raMulti1,1)
}
if ((raMulti1$bic[[2]][1]>1) && (raMulti1$bic[[2]][2])>1) {
message("\n\nPlot bicluster 2:
  1) bicluster in whole matrix
  2) only bicluster")
    plotBicluster(raMulti1,2)
}

message("\n\nPlot of pairs of biclusters as biplots:
  1) rectangles are samples
     color corresponds to tissue type
     desired: separation of the tissue types
  2) circles are genes
     red circles correspond to most indicative genes")
message("\n Plot1: Biclusters 1 and 2")
devAskNewPage(ask = TRUE)
plot(resMulti1,dim=c(1,2),label.tol=0.01,col.group=CMulti,lab.size=0.6)
message("\n Plot1: Biclusters 1 and 3")
plot(resMulti1,dim=c(1,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
message("\n Plot1: Biclusters 2 and 3")
plot(resMulti1,dim=c(2,3),label.tol=0.01,col.group=CMulti,lab.size=0.6)
devAskNewPage(ask = FALSE)
}
