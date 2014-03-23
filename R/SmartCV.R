##############################################################

#Proportional test subsetting for classification cv
PropSample <- function(ytemp,testprop=.25,nfold=10,minnum=3,output=F){
  ytemp=factor(ytemp)
  tbl=table(ytemp)
  badclass=unique(c(which(names(tbl)==""),which(names(tbl)=="Unassigned"),which(tbl<minnum)))
  names(badclass)=names(tbl)[badclass]
  if(sum(is.na(ytemp))>0){
    badind=which(is.na(ytemp))
  }else{
    badind=c()
  }
  if(length(badclass)>0){
    if(output==T){
      cat("PropSample Warning: Class(es)",paste(names(tbl)[badclass],collapse=",\n"),"has/have insufficient (<",minnum,") number of classes!\nOffending classes removed\n")
    }
    badnames=names(badclass)
    goodnames=names(tbl[-badclass])
    for(i in 1:length(badnames)){
      badind=c(badind,which(ytemp==badnames[i]))
    }
    goodind=c(1:length(ytemp))[-badind]
  }else{
    goodind=c(1:length(ytemp))
    goodnames=names(tbl)
  }
  
  #Checks to see if more trainind than nfold; if not, nfold lowered
  if(floor(length(goodind)*(1-testprop))<=nfold){
    nfold=floor(length(goodind)*(1-testprop))
  }
  #List holding test set indices and matrix for folds
  cvntest=list(NULL,matrix(0,nrow=nfold,ncol=max(1,ceiling(length(goodind)/nfold))))
  if(length(goodind)>0){
    matind=1
    ##Check ortho and soymo
    for(i in 1:length(goodnames)){
      #Collect indices of relevant label
      tind=which(ytemp==goodnames[i])
      #Randomize order of labels
      tind=tind[sample(1:length(tind),length(tind))]
      
      lasttestind=max(floor(length(tind)*testprop),1)
      testind=tind[1:lasttestind]
      trainind=tind[(lasttestind+1):length(tind)]
      cvntest[[1]]=c(cvntest[[1]],testind)
      cvntest[[2]][matind:(matind-1+length(trainind))]=trainind
      matind=matind+length(trainind)
      #     cat(matind,dim(cvntest[[2]]),"\n")
    }
  }
  cvntest2=vector("list",nfold+2)
  cvntest2[[1]]=cvntest[[1]]
  for(i in 1:nrow(cvntest[[2]])){
    cvntemp=cvntest[[2]][i,]
    cvntest2[[i+1]]=cvntemp[cvntemp!=0]
  }
  if(length(badind)>0){
    cvntest2[[nfold+2]]=badind
  }
  #First entry in the list is the test set, the rest are folds
  cvntest2
}

############################################################
# Transforms testing data with same transformations
# (boxcox, then scaling) as training data
TrainTestTrans=function(Train,Test,bc=T,scl=T,BCInstances=8,output=F){
  
  BCtrain=MVBoxCox(xmat=Train,lamseq=seq(-2,2,.05),NumInstances=BCInstances,Instance=1,output=output)
  Train=BCtrain$xt
  lambda=BCtrain$lambda
  smpv=BCtrain$smpv
  
  Train=scale(Train)
  center=attributes(Train)$'scaled:center'
  Scale=attributes(Train)$'scaled:scale'
  badcol=which(Scale==0)
  Train[,badcol]=0
  for(j in 1:length(lambda)){
    if(Scale[j]!=0){
      Test[,j]=Test[,j]+smpv[j]
      if(abs(lambda[j])<1e-10){
        Test[,j]=log(Test[,j])
      }else{
        Test[,j]=(Test[,j]^(lambda[j])-1)/lambda[j]
      }
    }
  }
  Test=t((t(Test)-center)/Scale)
  Test[,badcol]=0
  
  #####################################
  #Correcting outlying values
  trainmins=apply(Train,2,min)
  testmins=apply(Test,2,min)
  badmins=which(testmins<trainmins)
  if(length(badmins)>0){
    for(j in 1:length(badmins)){
      bmtmp=which(Test[,badmins[j]]<trainmins[badmins[j]])
      Test[bmtmp,badmins[j]]=trainmins[badmins[j]]
    }
  }
  
  trainmaxs=apply(Train,2,max)
  testmaxs=apply(Test,2,max)
  badmaxs=which(testmaxs>trainmaxs)
  if(length(badmaxs)>0){
    for(j in 1:length(badmaxs)){
      bmtmp=which(Test[,badmaxs[j]]>trainmaxs[badmaxs[j]])
      Test[bmtmp,badmaxs[j]]=trainmaxs[badmaxs[j]]
    }
  }    
  Test[is.na(Test)]=0
  
  
  list(Train=Train,Test=Test)
}
