# BiSET

The R version for developing BiSEAT toolbox is 4-1-3. Thus, please make sure that R version >= 4-1-3 when running BiSEAT.

In addition, please note that the "gWidgets2" package needs to be downloaded manually because it has been removed from the R language. Here, the "gWidgets2" package that needs to be downloaded manually has been uploaded to the folder "Packages". To avoid arising from unnecessary problem, please copy all the packages in the "Packages" folder to the R library before running the BiSEAT toolbox.

Install the package:

devtools::install_github('Chuheming/BiSET')

Start Running the package:


1,   library(BiSET)

2,   BiSET：：Frame()

BiSET can generate the synthetic dataset to verify the biclustering performance according to two evaluate method, Relevance and Recovery.

# Usage
#Generate the synthetic dataset

We consider that the biclustering partten, such as shift partten, scale partten and shift-partten, could be represented by one formular:

X = alpha + Beta * Gamma

Therefore, the synthetic dataset could be generated by setting three parameters and selecting the dataset whether has noise or overlap

#Biclustering
Seven classical biclustering algorihtms are collection in BiSEAT. They are CC, QUBIC, Bimax, rUnibic, BiBit, FABIA are Plaid. 

You could select the algorithm on the checkbox When you determin to use the biclustering algorithm to analysis your dataset. Then,
click the button 'Run'. Wait moment, a prompt box will alert the result.

Finally, you can click the button 'Save_Bic' to solve the result.

#GO Enrichment analysis and KEGG Analysis

This function will help you easily to get the enrichment result.
The format of the analyzed gene sample is as follows:

Symobol

A

B

C

...

Z

You only select the gene Symbol or ENSENBL and select the dataset specials, and also can customize the p-value. 
Then, click the button 'GO Analysis' or "KEGG Analysis".

Similarly, you also can save the enrichment resutl when you click the button 'Save_GO'.

