

Frame <- function(){
  #GUI
  library(tidyverse)
  library(stringr)
  library(readr)
  library(readxl)
  library(writexl)
  library(gWidgets2)
  library(gWidgets)
  library(gWidgetsRGtk2)
  #图
  library(dplyr)
  library(ggplot2)
  library(AnnotationDbi)
  #Biclustering algorithm
  library(biclust)
  library(BiBitR)
  library(runibic)
  library(QUBIC)
  library(fabia)
  #GO enrichment
  library(clusterProfiler)
  library(KEGG.db)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)

  real_path = all_data = data_source = NULL
  result_go = NULL
  # my_path_raw = NULL

  # th = gfile(text = "Select",type = "selectdir",initial.dir=getwd())


  #Ignore warning
  options(warn = -1)
  # Synthetic RGTK2 package and GTK Widget
  options(guiToolkit="RGtk2")
  #Synthetic main window
  win <- gwindow('R Window',visible = FALSE)

  #Synthetic toolbar
  # Synthetic open toolbar
  my_open <- gaction(label = "Open", icon = "open", handler = function(h,...){
    #Implement the handler function to open the GUI data that needs to be imported
    real_path = choose.files(caption = "Choose One File (.txt)")
    if(length(real_path)==0)
      galert('Please Select the Files to be Processed', title = "Tips", delay = 6)
    else{

      if (grepl("\\.txt$", real_path)){
        data_source_id <- as.matrix(read.table(real_path, sep="", dec=".", header=FALSE, stringsAsFactors=FALSE))
        galert('The data is loaded successfully!', title = "Tips",delay = 6)
      }else{
        galert('Please choose the correct file format', title = "Tips", delay = 6)
      }
    }
    data_source <<- data_source_id
  })

  my_about_content <- "The main function of BiSET is to customize the generation of simulated datasets. At the same time, the software provides six biclustering algorithms (CC,Bimax,QUBIC,FABIA,RUnibic,Plaid) to cluster the Syntheticd simulated datasets and real datasets. In addition, BiSET can also perform GO enrichment analysis."
  my_about <- gaction(label = "About", icon = 'about',handler = function(h,...){
    gmessage(my_about_content,title = "About BiSET",parent = win)
  })

  #the defined events are placed on the main window
  my_list <- list(
    open = my_open,
    seq = list(separator = TRUE),
    about = my_about,
    seq = list(separator = TRUE)
  )

  pl = 0

  gf1 <- gframe(horizontal = FALSE, container = win)
  bg_note1 <- (container = gf1)
  gl_note <- glabel("Generate the Simulation Dataset",container=bg_note1)
  bg_gl1 <- (container = gf1)
  ##### backgroud Matrix and bicluster ########
  bgm_Ol <- glabel("Backgroud matrix:",container = bg_gl1)
  bgmr_value <- gedit(text = "100", width = 6, container = bg_gl1)
  bgmc_value <- gedit(text = "100", width = 6, container = bg_gl1)


  bi_Ol <- glabel("Bicluster matrix:",container = bg_gl1)
  bir_value <- gedit(text = "10", width = 6, container = bg_gl1)
  bic_value <- gedit(text = "10", width = 6, container = bg_gl1)
  ###########      Overlapping and Noise   ######
  #Noise
  # N_Ol <- glabel("Noise_level:",container = bg_gl1)
  # N_value <- gedit(text = "0", width = 4, container = bg_gl1)

  N_bal <- glabel("Noise:",container = bg_gl1)
  N_value <- gWidgets::gdroplist(items = c('0','1'),selected = 1, container = bg_gl1)

  Pa_bal <- glabel("Partten:",container = bg_gl1)
  Pa_value <- gWidgets::gdroplist(items = c('Shift','Scale','Shift-Scale'),selected = 1, container = bg_gl1)



  #Overlapping
  gl_Ol <- glabel("Overlapping_level:",container = bg_gl1)
  Ol_value <- gWidgets::gdroplist(items = c('0','1'),selected = 1, container = bg_gl1)
  gseparator(horizontal = FALSE, container = bg_gl1)

  bg_gl2 <- (container = gf1)
  ### Dividing line
  gseparator(horizontal = TRUE,container = bg_gl2)
  gseparator(horizontal = TRUE,container = bg_gl2)

  gf2 <- gframe(horizontal = FALSE, container = win)
  bg_note1 <- (container = gf2)
  ## select contant
  bg_gl2 <- (container = gf2)
  bg_note2 <- (container = gf2)

  gl_note <- glabel("Biclustering Partten",container=bg_note1)
  #number of implanted bicluster
  num_m <- glabel("Num:",container = bg_gl2)
  num_v <- gedit(text = "10", width = 4, container = bg_gl2)
  #Alpha
  gl_alpha <- glabel("Alpha:",container = bg_gl2)
  alpha_value <- gedit(text = "0.1", width = 4, container = bg_gl2)
  gseparator(horizontal = FALSE, container = bg_gl2)
  #Beta
  gl_Beta <- glabel("Beta:",container = bg_gl2)
  Beta_value <- gedit(text = "0.1", width = 4, container = bg_gl2)
  gseparator(horizontal = FALSE, container = bg_gl2)
  #Gamma
  gl_Gamma <- glabel("Gamma:",container = bg_gl2)
  Gamma_value <- gedit(text = "0.1", width = 4, container = bg_gl2)
  gseparator(horizontal = FALSE, container = bg_gl2)
  ########################
  bg <- (container=gf1)



  gbutton("Generate", container=bg, handler=function(h,...){
    Pva = svalue(Pa_value)
    if (identical(Pva,'Shift')){
      Pva = 1
    }else{
      if (identical(Pva,'Scale')){
        Pva = 2
      }else{
        Pva = 3
      }
    }

    #Determine whether the build is successful
    bo_m = 0
    #Select the noise
    level_noise <- as.numeric(svalue(N_value))
    #Select the overlapp
    level_ov <- as.numeric(svalue(Ol_value))
    #get the value of background matrix and bilcuster
    size_bgmr <- as.numeric(svalue(bgmr_value))
    size_bgmc <- as.numeric(svalue(bgmc_value))
    size_bir <- as.numeric(svalue(bir_value))
    size_bic <- as.numeric(svalue(bic_value))


    #get the number of cluster
    numv = as.numeric(svalue(num_v))
    #####
    bgm <- matrix(c(0),nrow = size_bgmr, ncol = size_bgmc)
    bim <- matrix(c(0),nrow = size_bir, ncol = size_bic)


    av = as.numeric(svalue(alpha_value))
    bv = as.numeric(svalue(Beta_value))
    gv = as.numeric(svalue(Gamma_value))

    if (level_ov == 0){
      if (level_noise == 0){

        if (Pva==1){
          # Shift partten
          out = getshift(numv,av,bv,size_bir,size_bic,bim,bgm)
          bo_m = 1
          all_bim <<- out$bim

        } else{
          if (Pva == 2){
            # scale partten
            out = getscale(numv,gv,bv,size_bir,size_bic,bim,bgm)
            bo_m = 1
            all_bim <<- out$bim
          }else{
            #shift-scale partten
            out = getshiht_scale(numv,av,gv,bv,size_bir,size_bic,bim,bgm)
            bo_m = 1
            all_bim <<- out$bim
          }
        }

      }else{
        #has noise，non-overlapping
        for (i in 1:size_bgmr){
          bgm[i,] = rnorm(size_bgmc,mean = 0, sd = 1)
        }
        # implant bicluster
        if (Pva==1){
          # shift
          out = getshift(numv,av,bv,size_bir,size_bic,bim,bgm)
          bim <<- out$bim
          bgm <<- out$bgm
          all_bim <<- bim
          bo_m = 1
        } else{
          if (Pva == 2){
            # scale
            out = getscale(numv,gv,bv,size_bir,size_bic,bim,bgm)
            bim <<- out$bim
            bgm <<- out$bgm
            all_bim <<- bim
            bo_m = 1
          }else{
            #shift-scale partten
            out = getshiht_scale(numv,av,gv,bv,size_bir,size_bic,bim,bgm)
            bim <<- out$bim
            bgm <<- out$bgm
            all_bim <<- bim
            bo_m = 1
          }
        }
      }

    }else{

      if (level_noise == 0){
        #no noise ,has overlapping (only in shift-scale partten)
        # The overlapping interval defaults to two rows and five columns
        #Build clusters
        out = getover(gv,bv,av,size_bir,size_bic,bim,bgm,numv)
        bo_m = 1
        bim <<- out$bim
        bgm <<- out$bgm
        all_bim <<- bim

      }else{
        #noise , overlapping
        for (i in 1:size_bgmr){
          bgm[i,] = rnorm(size_bgmc,mean = 0, sd = 1)
        }
        out = getover(gv,bv,av,size_bir,size_bic,bim,bgm,numv)
        bim = out$bim
        bgm <<- out$bgm
        all_bim <<- bim
        bo_m <<- 1
      }

    }
    if (bo_m == 1){
      all_data <<- bgm
      all_bim <<- bim
      data_source <<- bgm
      pl <<- 1

      if (isEmpty(all_data)){
        galert('There is no output for the time',title = "File Save Failure",delay = 6)

      }else{

        # my_path_save = choose.dir(getwd(),"Choose a suitable folder")
        #my_path_save = paste0(th,"\\data_result.txt")
        #my_path_raw = paste0(th,"\\raw.txt")
        #write.table(all_data, file = my_path_save,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
        #write.table(bim, file = my_path_raw,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
        galert('Synthetic Data success! ',title = "File Save Success",delay = 6)
      }

   }


  })


  gbutton("Save_Cluster",container = bg_gl2,handler = function(h,...){
    th <<- gfile(text = "Select",type = "selectdir",initial.dir=getwd())
    if (pl==1){

      my_path_save = paste0(th,"\\data_result.txt")
      my_path_raw = paste0(th,"\\raw.txt")
      write.table(all_data, file = my_path_save,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
      write.table(all_bim, file = my_path_raw,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
      galert(paste0('Success! The results are preserved in ',my_path_save,'\\data_result.txt'),title = "File Save Success",delay = 6)

    }else{
      galert('There is no output for the time',title = "File Save Failure",delay = 6)
    }
  })


  ###
  gseparator(horizontal = TRUE,container = bg_gl2)
  gseparator(horizontal = TRUE,container = bg_gl2)

  gl_note <- glabel("Num: the number of implanted bicluster",container=bg_note1)
  gl_note <- glabel("x=Alpha+Beta*Gamma",container=bg_note1)



  gtb <- gWidgets::gtoolbar(toolbarlist = my_list, container = win)
  #gframe
  gf <- gframe(horizontal = FALSE, container = win)
  ## select contant
  bg_gl <- (container = gf)


  ## 一个筛选框
  gl_bal <- glabel("Biclustering Algorithm:",container = bg_gl)
  bal_value <- gWidgets::gdroplist(items = c('CC','Bimax','Qubic','rUnibic','Plaid','FABIA'),selected = 1, container = bg_gl)
  gseparator(horizontal = FALSE,container = bg_gl)

  #Button
  sbo <- 0 #
  sbf <- 0 #FABIA

  gbutton("Run",container = bg_gl,handler = function(h,...){

    loma <<- data_source
    if (identical(svalue(bal_value),'CC')){


      if (length(loma)==0){
        galert("Please add or synthetic dataset first")

      }else
      {
        gmessage("Wait")
        sbo <<- 1
        re <<- biclust(loma, method=BCCC(), delta=0.5, alpha=1, number=100)

        tt = re@Number
        galert(paste0("The bicluster number is ",tt))

      }
    }
    if(identical(svalue(bal_value),'Bimax')){


      if (length(loma)==0){
        galert("Please add or synthetic dataset first")
      }else{
        gmessage("Wait")
        sbo <<- 1
        biloma = binarize(loma)
        re <<- biclust(x=biloma,method=BCBimax(),minc=2,minr=2,number=10)
        tt = re@Number
        galert(paste0("The bicluster number is ",tt))
      }

    }
    if(identical(svalue(bal_value),'Qubic')){


      if (length(loma)==0){
        galert("Please add or synthetic dataset first")
      }else{
        gmessage("Wait")
        sbo <<- 1
        re <<- biclust::biclust(loma,method = BCQU(),r=1,q=0.06,c=0.95,o=100,f=1)
        tt = re@Number
        galert(paste0("The bicluster number is ",tt))

      }

    }
    if(identical(svalue(bal_value),'rUnibic')){

      if (length(loma)==0){
        galert("Please add or synthetic dataset first")
      }else{
        gmessage("Wait")
        sbo <<- 1
        re <<-  runibic(loma)
        tt = re@Number
        galert(paste0("The bicluster number is ",tt))

      }

    }
    if(identical(svalue(bal_value),'Plaid')){

      if (length(loma)==0){
        galert("Please add or synthetic dataset first")
      }else{
        gmessage("Wait")
        sbo <<- 1
        re <<- biclust(loma, method = BCPlaid())#plaid
        tt = re@Number
        galert(paste0("The bicluster number is ",tt))
      }

    }
    if(identical(svalue(bal_value),'FABIA')){

      if (length(loma)==0){
        galert("Please add or synthetic dataset first")
      }else
      {
        sbf <<- 1
        gmessage("Wait")
        fre <<- fabia(loma,100,0.1,1000)
        rb = extractBic(fre)
        lrb = length(rb)
        galert(paste0("The bicluster number is ",lrb))

      }

    }

  })

  gl_bal <- glabel("Input the index:",container = bg_gl)
  se_value <- gedit(text = "1", width = 4, container = bg_gl)
  gbutton("Save_Bic",container = bg_gl,handler = function(h,...){

    th = gfile(text = "Select",type = "selectdir",initial.dir=getwd())
    my_path_save = paste0(th,"\\Bicrow.txt")
    my_path_raw = paste0(th,"\\Biccol.txt")
    #ouput birow，bicol
    k = as.numeric(svalue(se_value))
    if(identical(svalue(bal_value),'FABIA')){
      rb = extractBic(fre)
      lrb = length(rb)

        bicluster = rb$bic[k,]
        asd = bicluster$biypn
        r = bicluster$bixn
        birow = matrix(c(0),nrow = 1, ncol = lrb)

        r = as.matrix(t(r))
        asd = as.matrix(t(asd))

        rd = dim(r)
        ra = dim(asd)

        for (i in 1:rd[2]){
          r[1,i] = gsub('[gene]','',r[1,i])
        }
        for (j in 1:ra[2]){
          asd[1,j] = gsub('[V]','',asd[1,j])
        }
        b = which(r==0)
        r = as.matrix(r[-b])
        r = t(r)

        b = which(asd==0)
        asd = as.matrix(asd[-b])
        asd = t(asd)

        write.table(r, file = my_path_save,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
        write.table(asd, file = my_path_raw,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
        galert(paste0('Success! The results are preserved in ',my_path_save),title = "File Save Success",delay = 6)


    }else{
      tt = re@Number
      num = dim(re@NumberxCol)
      rnum = dim(re@RowxNumber)

      d = matrix(c(0),nrow = 1,ncol = rnum[1])
      jj = 1
      for (i in 1:rnum[1]){
        if (identical(re@RowxNumber[i,k],TRUE)){
          d[1,jj] = i # gene location
          jj = jj + 1
        }
      }
      cc = matrix(c(0), nrow = 1, ncol = num[2])
      jj = 1

      for (i in 1:num[2]){
        if (identical(re@NumberxCol[k,i],TRUE)){
          cc[1,jj] = i; # condition location
          jj = jj+1
        }
      }
      b = which(d==0)
      d = as.matrix(d[-b])
      d = t(d)
      print(d)

      b = which(cc==0)
      cc = as.matrix(cc[-b])
      cc = t(cc)
      print(cc)
      write.table(d, file = my_path_save,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
      write.table(cc, file = my_path_raw,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
      galert(paste0('Success! The results are preserved in ',my_path_save),title = "File Save Success",delay = 6)

    }
  })


  ##
  bg_note <- (container = gf)
  gseparator(horizontal = TRUE,container = gf)
  gseparator(horizontal = TRUE,container = gf)

  #synthetic lable
  gl_note <- glabel("Verify the simulate dataset",container=bg_note)
  gbutton("Recovery",container = bg_note,handler = function(h,...){
    my_path_raw = paste0(th,"\\raw.txt")
    raw = as.matrix(read.table(my_path_raw, sep="", dec=".", header=FALSE, stringsAsFactors=FALSE))
    craw <<- raw
    impr <<- dim(raw)
    if (sbf == 1){
      if (identical(svalue(bal_value),'FABIA')){
        rb = extractBic(fre)
        lrb = length(rb)
        score1 = matrix(c(0),nrow = 1, ncol = lrb)
        for (i in 1:lrb){
          bicluster = rb$bic[i,]
          d = Get_data(bicluster,loma)

          cd = dim(d)
          if (is.null(cd)){
            score1[1,i] = 0
          }else{
            score = Get_sc(cd,impr,craw,d)

            score1[1,i] = score
          }
        }
        score = mean(score1)
        galert(paste0("The Recovery score is: ",score))

      }
    }else{
      if (sbo == 1){
        tt = re@Number
        if(tt>0){
          rownum = re@RowxNumber
          score1 = matrix(c(0),nrow = 1, ncol = tt)
          num = dim(re@NumberxCol)

          rnum = dim(re@RowxNumber)
          rownum = re@RowxNumber
          score1 = matrix(c(0),nrow = 1, ncol = tt)
          num = dim(re@NumberxCol)

          rnum = dim(re@RowxNumber)
          for (pn in 1:tt){
            d = loma[re@RowxNumber[,pn],re@NumberxCol[pn,]]
            d = as.matrix(re@RowxNumber[,1])
            #
            cd = dim(d)
            score = Get_sc(cd,impr,craw,d)
            score1[1,pn] = score
          }
          score = mean(score1)
          galert(paste0("The Recovery score is: ",score))
        }else{
          galert("No Bicluster")
        }

      }else
        galert('No found file')
    }


  })
  gbutton("Relevance",container = bg_note,handler = function(h,...){
    my_path_raw = paste0(th,"\\raw.txt")
    raw = as.matrix(read.table(my_path_raw, sep="", dec=".", header=FALSE, stringsAsFactors=FALSE))
    craw <<- raw
    impr <<- dim(raw)
     if (sbf == 1){
      if (identical(svalue(bal_value),'FABIA')){
        rb = extractBic(fre)
        lrb = length(rb)
        score1 = matrix(c(0),nrow = 1, ncol = lrb)
        for (i in 1:lrb){
          bicluster = rb$bic[i,]
          d = Get_data(bicluster,loma)
          cd = dim(d)
          if (is.null(cd)){
            score1[1,i] = 0
          }else{
            score = Get_sc(cd,impr,d,craw)

            score1[1,i] = score
          }

        }
        score = max(score1)
        galert(paste0("The Relevance score is: ",score))

      }
    }else{

      if (sbo == 1){
        tt = re@Number
        if (tt>0){
          rownum = re@RowxNumber
          score1 = matrix(c(0),nrow = 1, ncol = tt)
          num = dim(re@NumberxCol)

          rnum = dim(re@RowxNumber)
          for (pn in 1:tt){
            d = loma[re@RowxNumber[,pn],re@NumberxCol[pn,]]

            d = as.matrix(re@RowxNumber[,1])

            cd = dim(d)
            score = Get_sc(cd,impr,d,craw)

            score1[1,pn] = score
          }
          score = max(score1)
          galert(paste0("The Relevance score is: ",score))
        }else{
          galert("No bicluster")
        }

      }else{
        galert('No found data')
      }
    }

  })

  ##################################################
  ### Go enrichment
  ###
  ###
  real_note <- (container = gf)
  gseparator(horizontal = TRUE,container = real_note)


  gl_bal <- glabel("GO Enrichment Analysis:",container = real_note)
  gl_bal <- glabel("Please select:",container = real_note)
  sy_value <- gWidgets::gdroplist(items = c('Symbol','ENSEMBL'),selected = 1,  container = real_note)
  gl_bal <- glabel("Please select:",container = real_note)
  spe_value <- gWidgets::gdroplist(items = c('Homo','Mus'),selected = 1,  container = real_note)
  gl_bal <- glabel("Please select:",container = real_note)
  p_value <- gedit(text = "0.05", width = 4, container = real_note)

  gbutton("GO Analysis",container = real_note,handler = function(h,...){

    syva = svalue(sy_value)
    if (identical(syva,'Symbol')){
      syva = 1
    }else{
      syva = 0
    }
    speva = svalue(spe_value)
    if (identical(speva,'Homo')){
      speva = 1
    }else
      speva = 0

    pvalue = as.numeric(svalue(p_value))
    acount = 0 #calculate the number of gross enrich cluster
    gene_path = choose.files(caption = "Choose One File (.txt)")
    mat = read.table(gene_path, sep="", dec=".", header=FALSE, stringsAsFactors=FALSE)

    dimmat = dim(mat)

    genesymbol = matrix(c(0),nrow = dimmat[1],ncol = 1)
    count = 0
    for (i in 1:dimmat[1]){
      if (identical(mat[i,1],"Symbol")){
        count = count + 1
      }
    }
    kcount = count
    location_sybom = matrix(c(0),nrow = kcount,ncol = 1)
    k = 1
    for (i in 1:dimmat[1]){
      if (identical(mat[i,1],"Symbol")){
        location_sybom[k,1] = i
        k = k + 1
      }
    }
    sizedimmat = dim(location_sybom)
    total = sizedimmat[1]

    k = 1
    facount = 1
    gsl = 1
    i = 1
    jjjj = 1
    while (i < dimmat[1]){

      if (sizedimmat[1]==1){
        genesymbol = matrix(c(0),nrow = dimmat-1,ncol = 1)

        for (jjj in 2:dimmat[1]){
          genesymbol[jjjj,1] = mat[jjj,1]
          jjjj = jjjj + 1
        }

       i = i+ jjjj

       result_go1 <- executGO(genesymbol,speva,syva,pvalue)
       derego1 <- dim(result_go1)
       if (derego1[1]>0){
         d_cbind <- rbind(result_go,result_go1)
         result_go <<- d_cbind
         acount = acount + 1
       }


      }else{
        ti = i+1
        kt = k+1
        if (kt<=total){
          sizek = location_sybom[kt,1]-location_sybom[k,1]-1
          genesymbol = matrix(c(0),nrow = sizek,ncol = 1)
          gsl = 1
          while (ti<location_sybom[kt,1]){
            genesymbol[gsl,1] = mat[ti,1]
            gsl = gsl + 1
            ti = ti + 1
          }
        }
        else{
          sizek = dimmat[1]-location_sybom[k,1]+1
          genesymbol = matrix(c(0),nrow = sizek,ncol = 1)
          gsl = 1
          while (ti<dimmat[1]){
            genesymbol[gsl,1] = mat[ti,1]
            gsl = gsl + 1
            ti = ti + 1
          }
        }
        k = k + 1
        i = ti
         result_go1 <<- executGO(genesymbol,speva,syva,pvalue)
         derego1 <- dim(result_go1)
         if (derego1[1]>0){
           d_cbind <- rbind(result_go,result_go1)
           result_go <<- d_cbind
           acount = acount + 1
         }

      }

    }

    result_go <<- result_go
    galert(paste0("Total: ",total,";","Enrich: ",count,"."))

  })

  gbutton("KEGG Analysis",container = real_note,handler = function(h,...){
    syva = svalue(sy_value)
    if (identical(syva,'Symbol')){
      syva = 1
    }else{
      syva = 0
    }
    speva = svalue(spe_value)
    if (identical(speva,'Homo')){
      speva = 1
    }else
      speva = 0

    pvalue = as.numeric(svalue(p_value))

    gene_path = choose.files(caption = "Choose One File (.txt)")
    mat = as.matrix(read.table(gene_path, sep="", dec=".", header=FALSE, stringsAsFactors=FALSE))

    gene <- unique(mat)

    if (speva == 1){
      if (syva == 1){
        gene.df <- bitr(gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
      }else
        gene.df <- bitr(gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)
      gene = gene.df$ENTREZID
      kk <- enrichKEGG(gene = gene, organism = "hsa",pvalueCutoff = 1, qvalueCutoff = 1,use_internal_data = TRUE)

    }else{
      if (syva == 1){
        gene.df <- bitr(gene,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)
      }else{
        gene.df <- bitr(gene,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Mm.eg.db)

      }
      gene = gene.df$ENTREZID
     kk <- enrichKEGG(gene = gene_up, organism = "mmu",pvalueCutoff = 1, qvalueCutoff = 1,use_internal_data = TRUE)

    }

    pic <- dotplot(kk)

    plot(pic)

  })
  gbutton("Save_GO",container = real_note,handler = function(h,...){
    if (!is.null(result_go)){
      th = gfile(text = "Select",type = "selectdir",initial.dir=getwd())
      my_path_save = paste0(th,"\\GOENRICH.txt")
      write.table(result_go, file = my_path_save,append = FALSE, quote = TRUE, sep = " ",eol = "\n", na = "NA", dec = ".", row.names = F, col.names = F)
      galert("Output Complete")
    }else
      galert("No Data output")

  })
  #################################################
  ##Set the size of the main window
  size(win) <- c(700, 650)
  ##main window visible
  visible(win) <- TRUE


}

Get_sc = function(cd,impr,d1,raw1){
  #
  cd = dim(d1)

  Sc = matrix(c(0),nrow=1,ncol=cd[1])
  for (i in 1 : cd[1]){
    sc = 0
    unsc = 0
    for (j in 1:impr[1])
    {
      ints = intersect(d1[i,],raw1[j,])
      s = dim(as.matrix(ints))
      intsr = s[2]

      if ( sc < intsr){
        sc = intsr
        nsc = dim(as.matrix(union(d1[i,],raw1[j,])))

        unsc = nsc[2]
      }
    }
    if (unsc>0){
      Sc[1,i] = sc/unsc

    }
  }
  X = sum(Sc)/cd[1]
  return (X)
}

Get_data = function(bicluster,loma){

  asd = bicluster$biypn
  r = bicluster$bixn
  if (!is_empty(asd) && !is_empty(r)){
    r = as.matrix(t(r))
    asd = as.matrix(t(asd))

    rd = dim(r)
    ra = dim(asd)

    for (i in 1:rd[2]){
      r[1,i] = gsub('[gene]','',r[1,i])
    }
    for (j in 1:ra[2]){
      asd[1,j] = gsub('[V]','',asd[1,j])
    }
    r = as.matrix(t(as.numeric(r)))
    asd = as.matrix(t(as.numeric(asd)))
    Data = matrix(c(0),nrow = rd[2],ncol = ra[2])
    for (i in 1:rd[2]){
      for (j in 1:ra[2]){
        Data[i,j] = loma[r[1,i],asd[1,j]]
      }
    }
  }else{
    Data = 0
  }
  return (Data)
}


executGO = function(genesymbol,speva,syva,pvalue){
  pvalueFilter = 0.05
  if (speva == 1){
    if (syva == 1){
      gene.df <- bitr(genesymbol,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求
    }else
      gene.df <- bitr(genesymbol,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Hs.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求
    gene = gene.df$ENTREZID
    kk = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", readable = T)
    Go_cc = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="CC", readable = T)
    Go_bp = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="BP", readable = T)
    Go_mf = enrichGO(gene = gene, OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="MF", readable = T)

  }else{
    if (syva == 1){
      gene.df <- bitr(genesymbol,fromType="SYMBOL",toType="ENTREZID", OrgDb = org.Mm.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求
    }else{
      gene.df <- bitr(genesymbol,fromType="ENSEMBL",toType="ENTREZID", OrgDb = org.Mm.eg.db)#TCGA数据框如果没有进行基因注释，那么fromType应该是Ensembl，各种ID之间可以互相转换,toType可以是一个字符串，也可以是一个向量，看自己需求

    }
    gene = gene.df$ENTREZID
    kk = enrichGO(gene = gene, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont = "all", readable = T)
    Go_cc = enrichGO(gene = gene, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="CC", readable = T)
    Go_bp = enrichGO(gene = gene, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="BP", readable = T)
    Go_mf = enrichGO(gene = gene, OrgDb = org.Mm.eg.db, pvalueCutoff = 0.05, qvalueCutoff = 0.05, ont="MF", readable = T)

  }
  GO_all_result = as.data.frame(kk)
  GO_result_cc = Go_cc[(Go_cc$pvalue<pvalueFilter),]
  cc_num = dim(GO_result_cc)
  cc_num = cc_num[1]
  GO_result_mf = Go_mf[(Go_mf$pvalue<pvalueFilter),]
  mf_num = dim(GO_result_mf)
  mf_num = mf_num[1]
  GO_result_bp = Go_bp[(Go_bp$pvalue<pvalueFilter),]
  bp_num = dim(GO_result_bp)
  bp_num = bp_num[1]
  Go_all_result_enrich=GO_all_result[(GO_all_result$pvalue<pvalueFilter),]
  all_num = dim(Go_all_result_enrich)
  all_num = all_num[1]
  if(!is.null(Go_all_result_enrich)){
    galert("Drawing, please wait",delay = 6)
    ##
    display_number = c(20,20,20)
    BP_result <- as.data.frame(Go_bp)[1:display_number[1],]
    CC_result <- as.data.frame(Go_cc)[1:display_number[2],]
    MF_result <- as.data.frame(Go_mf)[1:display_number[3],]

    go_enrich_df <- data.frame(ID = c(BP_result$ID,CC_result$ID,MF_result$ID),
                               Description=c(BP_result$Description,CC_result$Description,MF_result$Description),
                               GeneNumber = c(BP_result$Count,CC_result$Count,MF_result$Count),
                               type=factor(c(rep("BP",display_number[1]),rep("CC",display_number[2]),rep("MF",display_number[3])),levels = c("BP","CC","MF")))

    for (i in 1:nrow(go_enrich_df)){
      description_splite = strsplit(go_enrich_df$Description[i],split = " ")
      description_collapse = paste(description_splite[[1]][1:5],collapse = " ")
      go_enrich_df$Description[i] = description_collapse
      go_enrich_df$Description = gsub(pattern = "NA","",go_enrich_df$Description)
    }

    #
    go_enrich_df$type_order = factor(rev(as.integer(rownames(go_enrich_df))), labels=rev(go_enrich_df$Description))
    COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")


    gg1 <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) + #Horizontal and vertical axis values
                       geom_bar(stat="identity", width=0.8) + #the width of the histogram
                       scale_fill_manual(values = COLS) + ###color
                       coord_flip() + ##Make the bar chart sideways
                       xlab("GO term") +
                       ylab("Gene_Number") +
                       labs(title = "The Most Enriched GO Terms")+
                       theme_bw()

    print(gg1)
    go_enrich_df$type_order=factor(rev(as.integer(rownames(go_enrich_df))),labels=rev(go_enrich_df$Description))

    COLS <- c("#66C3A5", "#8DA1CB", "#FD8D62")
    gg <- ggplot(data=go_enrich_df, aes(x=type_order,y=GeneNumber, fill=type)) +
                      geom_bar(stat="identity", width=0.8) +
                      scale_fill_manual(values = COLS) +
                      theme_bw() +
                      xlab("GO term") +
                      ylab("Num of Genes") +
                      labs(title = "The Most Enriched GO Terms")+
                      theme(axis.text.x=element_text(face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 ))

    print(gg)
    gmessage(paste0("GO enrichment items: ",all_num,";","CC items: ",cc_num,";","BP items: ",bp_num,";","MF items: ",mf_num,"."),title = "GO enrichment Information")
  }else{
    galert("No enrich item")
  }
  return(Go_all_result_enrich)
}

getshift = function(numv,av,bv,size_bir,size_bic,bim,bgm){
  bo_m= 1

  for (i in 1:numv){

    av1 = av
    m = 1
    m = m+(i-1)*size_bir
    mr = m+size_bir-1
    mc = m+size_bic-1
    for (k in m:mr){
      bv1 = bv
      for (kk in m:mc){
        bgm[kk,k] = av1 + bv1
        bv1 = bv1 + bv1*0.2
      }
      av1 = av1 + av1*0.2
    }
  }

  #build cluster

  av1 = av
  for(i in 1:size_bir){
    bv1 = bv
    for (j in 1:size_bic){
      bim[j,i] = av1 + bv1
      bv1 = bv1 + bv1*0.2
    }
    av1 = av1 + av1*0.2
  }
  return(list(bim = bim,bgm = bgm))
}

getscale = function(numv,gv,bv,size_bir,size_bic,bim,bgm){
  bo_m = 1
  for (i in 1:numv){
    m = 1
    m = m+(i-1)*size_bir
    mr = m+size_bir-1
    mc = m+size_bic-1
    gv1 = gv

    for (k in m:mr){
      bv1 = bv
      for (kk in m:mc){
        bgm[kk,k] = bv1 * gv1
        gv1 = gv1 + gv1*0.02
      }
      bv1 = bv1 + bv1*0.02
    }
  }
  #Build clusters
  gv1 = gv

  for(i in 1:size_bir){
    bv1 = bv
    for (j in 1:size_bic){
      bim[j,i] = bv1 * gv1
      gv1 = gv1 + gv1*0.2
    }
    bv1 = bv1 + bv1*0.1
  }
  return(list(bim = bim,bgm = bgm))
}

getshiht_scale = function(numv,av,gv,bv,size_bir,size_bic,bim,bgm){
  bo_m = 1
  for (i in 1:numv){
    m = 1
    m = m+(i-1)*size_bir
    mr = m+size_bir-1
    mc = m+size_bic-1

    gv1 = gv
    for (k in m:mr){
      bv1 = bv
      av1 = av
      for (kk in m:mc){
        bgm[kk,k] = av1 + bv1 * gv1
        av1 = av1 + av1 * 0.05
        bv1 = bv1 + bv1 * 0.2
      }
      gv1 = gv1 + gv1*0.02
    }
  }
  #Build clusters

  gv1 = gv

  for(i in 1:size_bir){
    bv1 = bv
    av1 = av
    for (j in 1:size_bic){
      bim[j,i] = av1 + bv1 * gv1
      av1 = av1 + av1 * 0.05
      bv1 = bv1 + bv1 * 0.2
    }
    gv1 = gv1 + gv1*0.02
  }
  return(list(bim = bim,bgm = bgm))
}

getover = function(gv,bv,av,size_bir,size_bic,bim,bgm,numv){
  gv1 = gv

  for(i in 1:size_bir){
    bv1 = bv
    av1 = av
    for (j in 1:size_bic){
      bim[j,i] = av1 + bv1 * gv1
      av1 = av1 + av1 * 0.05
      bv1 = bv1 + bv1 * 0.2
    }
    gv1 = gv1 + gv1*0.02
  }
  all_bim <<- bim

  fbim = matrix(c(0),nrow = size_bir,ncol = size_bic)
  for (i in 1:size_bir){
    temp_i = size_bir-i+1
    fbim[i,] <- bim[temp_i,]
  }

  bo_m = 1
  for (i in 1:numv){
    m = 1

    #column start position
    cm = m+(i-1)*(size_bic-5)
    #line start position
    m = m+(i-1)*(size_bir-2)
    #line end positon
    mr = m+size_bir-1
    fi1 = 3
    for (k in 1:size_bir){
      if (i%%2==0){
        #odd i put fbim

        if (k<=2){
          jc <- cm+5
          jcc = cm + size_bic-1
          kk = 2-k+1
          bgm[m,jc:jcc] = fbim[kk,1:5]
          m = m + 1
        }else{
          j2 = cm + size_bir-1
          bgm[m,cm:j2] = fbim[fi1,1:size_bic]
          m = m + 1
          fi1 = fi1 + 1
        }


      }else{
        #
        if (i>1){

          if (k<=2){
            jc <- cm+5
            jcc = cm + size_bic-1
            kk = 2-k+1
            sc = size_bic-5
            bgm[m,jc:jcc] = bim[kk,1:sc]
            m = m + 1
          }else{
            j2 = cm + size_bir-1
            bgm[m,cm:j2] = bim[fi1,1:size_bic]
            m = m + 1
            fi1 = fi1 + 1
          }

        }else{
          trj2 = cm+ size_bic-1

          bgm[k,cm:trj2] = bim[k,]
        }
      }
    }
  }
  return(list(bim = bim,bgm = bgm))
}
