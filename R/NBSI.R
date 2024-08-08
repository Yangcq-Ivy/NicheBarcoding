

#' Niche-model-Based Species Identification (NBSI)
#'
#' @description Species identification using DNA barcoding integrated with
#' niche model.
#'
#' @param  ref.seq DNAbin, the reference dataset containing sample IDs,
#' taxon information,longitude and latitude, and barcode sequences of samples.
#' @param  que.seq DNAbin, the query dataset containing sample IDs, longitude
#' and latitude, and barcode sequences of samples.
#' @param  model Character, string indicating which niche model will be used.
#' Must be one of "MAXENT" (default) or "RF". "MAXENT" can only be applied when
#' the java program paste(system.file(package="dismo"), "/java/maxent.jar",
#' sep='') exists.
#' @param  ref.add Data.frame, the additional coordinates collected from GBIF
#' or literatures.
#' @param  en.vir RasterBrick, the global bioclimate data output from
#' "raster::getData" function.
#' @param  bak.vir Matrix, bioclimate variables of random background points.
#'
#' @return A dataframe of barcoding identification result for each query sample
#' and corresponding niche model-based probability.
#'
#' @keywords NBSI
#' @export
#'
#' @author Cai-qing YANG (Email: yangcq(at)cnu.edu.cn) and Ai-bing ZHANG
#' (Email:zhangab2008(at)cnu.edu.cn), Capital Normal University (CNU), Beijing,
#' CHINA.
#'
#' @references Breiman, L. 2001. Random forests. Machine Learning 45(1):5-32.
#' @references Liaw, A. and M. Wiener. 2002. Clasification and regression by
#' randomForest. R News, 2/3:18-22.
#' @references Phillips, S.J., R.P. Anderson and R.E. Schapire. 2006. Maximum
#' entropy modeling of species geographic distributions. Ecological Modelling,
#' 190:231-259.
#' @references Zhang, A.B., M.D. Hao, C.Q. Yang and Z.Y. Shi. (2017). BarcodingR:
#' an integrated R package for species identification using DNA barcodes.
#' Methods in Ecology and Evolution, 8:627-634.
#' @references Jin, Q., H.L. Han, X.M. Hu, X.H. Li, C.D. Zhu, S.Y.W. Ho, R.D. Ward
#' and A.B. Zhang. 2013. Quantifying species diversity with a DNA barcoding-based
#' method: Tibetan moth species (Noctuidae) on the Qinghai-Tibetan Plateau.
#' PloS One, 8:e644.
#' @references Hijmans, R.J., S.E. Cameron, J.L. Parra, P.G. Jones and A. Jarvis.
#' 2005. Very high resolution interpolated climate surfaces for global land areas.
#' International Journal of Climatology, 25(15):1965-1978.
#'
#'
#' @examples
#' data(en.vir)
#' data(bak.vir)
#' #envir<-raster::getData("worldclim",download=FALSE,var="bio",res=2.5)
#' #en.vir<-raster::brick(envir)
#' #back<-dismo::randomPoints(mask=en.vir,n=5000,ext=NULL,extf=1.1,
#' #                          excludep=TRUE,prob=FALSE,
#' #                          cellnumbers=FALSE,tryf=3,warn=2,
#' #                          lonlatCorrection=TRUE)
#' #bak.vir<-raster::extract(en.vir,back)
#'
#' library(ape)
#' data(LappetMoths)
#' ref.seq<-LappetMoths$ref.seq[1:50,]
#' que.seq<-LappetMoths$que.seq[1:5,]
#' NBSI.out<-NBSI(ref.seq,que.seq,ref.add=NULL,
#'                model="RF",en.vir=en.vir,bak.vir=bak.vir)
#' NBSI.out
#'
#' ### Add a parameter when additional reference coordinates are available ###
#' #ref.add<-LappetMoths$ref.add
#' #NBSI.out2<-NBSI(ref.seq,que.seq,ref.add=ref.add,
#' #                model="RF",en.vir=en.vir,bak.vir=bak.vir)
#' #NBSI.out2


NBSI<-function(ref.seq,que.seq,model="MAXENT",ref.add=NULL,
               en.vir=NULL,bak.vir=NULL)
{
  ### 1 Ref.mono.test ####
  ref.infor<-extractSpeInfo(rownames(ref.seq))
  que.infor<-extractSpeInfo(rownames(que.seq))
  que.var <- raster::extract(en.vir,que.infor[,4:5])
  que.var <- as.data.frame(que.var)
  if (all(is.na(que.var)) == TRUE){
    stop("No variables can be extracted from que.infor!")
  }
  row.names(ref.seq)<-gsub("\\,[0-9\\.\\ \\-]*$","",rownames(ref.seq));ref.seq
  row.names(que.seq)<-gsub("\\,[0-9\\.\\ \\-]*$","",rownames(que.seq));que.seq

  sampleSpeNames <- attr(ref.seq, "dimnames")[[1]]
  Spp <- gsub(".+,", "", sampleSpeNames)
  refseq.tree <- ape::nj(stats::dist(ref.seq))
  seq.mono <- monophyly.prop(refseq.tree, Spp, singletonsMono = FALSE)
  mono.list.refseq <- seq.mono$mono.list

  ref.env <- raster::extract(en.vir,ref.infor[,4:5])
  rownames(ref.env) <- ref.infor[,3]
  refenv.tree <- ape::nj(stats::dist(ref.env))
  env.mono <- monophyly.prop(refenv.tree, ref.infor[,3], singletonsMono = FALSE)
  mono.list.refenv <- env.mono$mono.list

  ### 2 Bayesian species identification ####
  Bayesian.all.prob<-function(ref,que,all.ref.sp){
    ref2<-as.data.frame(as.character(ref))
    que2<-as.data.frame(as.character(que))

    ref3<-cbind(ref2,species=gsub(".+,","",rownames(ref)))
    rownames(ref3)<-1:nrow(ref3)
    que3<-cbind(que2,species=gsub(".+,","",rownames(que)))

    Bayesian.trained <- e1071::naiveBayes(species ~ ., data = ref3)
    spe.inferred<-stats::predict(Bayesian.trained, que3)
    spe.inferred.prob<-stats::predict(Bayesian.trained, que3, type = "raw")

    all.prob<-t(spe.inferred.prob)
    all_prob<-list()
    for (sip in 1:nrow(que)){
      all_prob[[sip]]<-data.frame(usp=rownames(all.prob),tmp.sip=all.prob[,sip])
      rownames(all_prob[[sip]])<-1:length(all.ref.sp)
    }
    all_prob

    Bayesian.prob<-apply(spe.inferred.prob,1,max)

    spe.inferred<-as.character(spe.inferred)

    output.identified<-cbind(queryID=rownames(que),
                             species.identified=as.character(spe.inferred),
                             barcoding.based.prob=Bayesian.prob)
    output.identified<-as.data.frame(output.identified)

    out<-list(output_identified=output.identified,all_prob=all_prob)


    class(out) <- c("BarcodingR")
    return(out)
  }

  all.ref.sp<-unique(ref.infor[,3])
  bsi0<-Bayesian.all.prob(ref.seq,ref.seq,all.ref.sp)
  ref.res <- bsi0$output_identified
  TP=0;TN=0;FP=0;FN=0
  for (rrs in 1:nrow(ref.res)){
    q <- gsub(".+,","",as.character(ref.res[rrs,1]))
    out <- as.character(ref.res[rrs,2])
    fmf <- ref.res[rrs,3]
    if (q == out){
      ifelse(fmf >= 0.99, (TP=TP+1), (FN=FN+1))
    }else{
      ifelse(fmf >= 0.99, (FP=FP+1), (TN=TN+1))
    }
  }
  success.b=(TP+TN)/nrow(ref.res);success.b
  barcode.w <- success.b

  bsi<-Bayesian.all.prob(ref.seq,que.seq,all.ref.sp)
  BSI.out<-bsi$output_identified
  BSI.out[,1]<-gsub(",unknown","",BSI.out[,1])
  BSI.out$barcoding.based.prob<-as.numeric(BSI.out$barcoding.based.prob)
  all.ref.sp <- unique(BSI.out[,2])

  ### 3 Niche Modeling ####
  potent.ref.infor <- ref.infor[ref.infor[,3] %in% all.ref.sp,]
  samp.env <- raster::extract(en.vir,potent.ref.infor[,4:5])
  species <- potent.ref.infor[,3]

  eff.samp.env<-as.data.frame(cbind(species,samp.env))

  model<-gsub("randomforest|RandomForest|randomForest","RF",model)
  model<-gsub("maxent|Maxent","MAXENT",model)
  colnames(ref.infor)[3]<-gsub("Species|SPECIES","species",colnames(ref.infor)[3])

  # (1) Niche modeling and self-inquery of ref
  spe.niche <- list()
  sel.env <- list()
  niche.ref.w <- list()
  niche.ref.prob <- list()
  for (siu in 1:length(all.ref.sp)){
    prese.env<-eff.samp.env[gsub(".+,","",eff.samp.env[,1]) %in% all.ref.sp[siu],-1]
    prese.env<-stats::na.omit(prese.env)
    if (dim(prese.env)[1]==1){
      prese.env<-rbind(prese.env,prese.env)
      warning ("The model may not be accurate because there is only one record of ",all.ref.sp[siu],"!\n")
    }

    mod<-niche.Model.Build(prese=NULL,absen=NULL,
                           prese.env=prese.env,absen.env=NULL,
                           model=model,bak.vir=bak.vir,en.vir=en.vir)
    spe.niche[[siu]]<-mod$model
    sel.env[[siu]] <- mod$SelectedVariables

    spe.var<-apply(prese.env,FUN=as.numeric,MARGIN=2)
    spe.var<-as.data.frame(spe.var)

    success.n <- mod$SST[4] #Accuracy
    niche.ref.prob[[siu]]<-mod$SST[3] #Threshold
    niche.ref.w[[siu]] <- success.n
  }
  sel.env
  niche.ref.prob
  niche.ref.w

  # (2) Prediction of query samples and calculation of Prob(Sid)
  result <- data.frame()
  for(q in 1:nrow(que.seq)){
    queID<-as.character(que.infor[q,2])
    #target.spe<-as.character(BSI.out$output_identified[q,2]) # for FuzzyID
    target.spe<-as.character(BSI.out[q,2])  # for Bayesian
    if (length(target.spe) == 0){
      stop(paste("Please check the barcoding identification of ", queID,"!\n",sep=""))
      print(BSI.out)
    }
    spe_in_ref<-as.character(ref.infor[grep(target.spe,ref.infor$species),]$species)

    if (length(spe_in_ref) == 0){
      warning ("The identified species ",target.spe, " doesn't exist in ref.infor!",
               " Skipping the niche-based procedure ", queID," ...\n")

    }else{
      spe.index<-grep(paste(target.spe,"$",sep=""),all.ref.sp,fixed=F)
      model.spe<-spe.niche[[spe.index]]  #the model of target species

      if (all(is.na(que.infor[,4:5])) == TRUE){
        stop ("Please check the coordinate of ",queID," !\n")
      }else{
        que.var <- raster::extract(en.vir,que.infor[q,4:5])
        que.var <- as.data.frame(que.var)
        if (model == "RF"){
          que.HSI<-stats::predict(model.spe,que.var)
        }else if (model == "MAXENT"){
          que.HSI<-dismo::predict(model.spe,que.var,args='outputformat=logistic')
        }

        niche.w<-niche.ref.w[[spe.index]]
        niche.rocT<-niche.ref.prob[[spe.index]]

        # Find the nearest neighbor
        otherSP.env <- ref.env[-which(rownames(ref.env) == target.spe),sel.env[[spe.index]]]
        if (model == "RF"){
          otherSP.HSI<-stats::predict(model.spe,otherSP.env)
        }else if (model == "MAXENT"){
          otherSP.HSI<-dismo::predict(model.spe,otherSP.env,args='outputformat=logistic')
        }
        theta2 = 1-max(otherSP.HSI)
        theta1 = 1-niche.rocT
        x = 1-que.HSI

        # Calculate the Fuzzy probability of Pe
        if (x <= theta1){
          Pe = 1
        }else if (x >= theta1 && x <= (theta1+theta2)/2){
          Pe = 1-2*((x-theta1)/(theta2-theta1))^2
        }else if (x >= (theta1+theta2)/2 && x <= theta2){
          Pe = 2*((x-theta2)/(theta2-theta1))^2
        }else if (x >= theta2){
          Pe = 0
        }

        # Calculation of probability
        #Pb <- BSI.out$output_identified[q,3] # for FuzzyID
        Pb <- BSI.out[q,3]  # for Bayesian
        if (Pb >= 0.98 && Pe == 1){ #二者都同意
          NicoB.prob = 1
        }else if (Pb < 0.95 && Pe == 0){ #二者都不同意
          NicoB.prob = 0
        }else{ #二者出现分歧
          if (length(grep(BSI.out[q,2], mono.list.refseq)) != 0){ #如果seq单系性好
            w.sum = barcode.w + niche.w
            NicoB.prob = (barcode.w/w.sum)*Pb + (niche.w/w.sum)*Pe
          }else{  #如果seq单系性不好
            NicoB.prob = Pe
          }
        }
      }
      res0 <- cbind(Pb, Pe, NicoB.prob)
      res0<-cbind(queID,que.infor[q,3],target.spe,round(res0,4))
      result <- rbind(result, res0)
    }
  }
  colnames(result) <- c("queID", "prior.spe", "target.spe",
                        "Pb", "Pe", "NBSI.prob")
  rownames(result) <- c(1:nrow(que.seq))

  return(result)
}

# The end of NBSI #
