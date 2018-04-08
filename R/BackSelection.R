#' Build the model through backward selection
#'
#'
#' @param data A data frame in the long (interval) format with one line per unit of time.
#' @param Type A vector consisting the name of variables representing the start and stop of each
#'             time interval and event indicator. e.g. c("start","stop","event").
#' @param variables A vector consisting the name of varialbes that will be adjusted in the model.
#' @param continuous A vector of binary indicators representing whether the corresponding variables
#'                   are continuous variables (1=Yes, 0=NO)
#' @param TD A vector representing whether to force the corresponding variable having a time-dependent
#'           effect. Only takes three values:\cr
#'           1 = force TD effect of the corresponding variable\cr
#'           0 = do not force any effect to the corresponding variable \cr
#'           -1 = force the PH effect of the corresponding variable \cr
#' @param NL A vector representing whether to force the corresponding variable having a non-linear
#'           effect. Only takes three values:\cr
#'           1 = force NL effect of the corresponding variable \cr
#'           0 = do not force any effect to the corresponding variable \cr
#'           -1 = force the linear effect of the corresponding variable \cr
#' @param m The number of interior knots used in the regression B-spline, default value is 1
#' @param p The degree of the regression B-spline, default value is 2
#' @param alpha_back significane level used in the backward elemination process
#' @param knots Default value is -999
#' @note Note that the TD and NL arguments have different meanings than that in the \code{\link{CoxFlex}} and \code{\link{tvcFlex}}
#' @export




BackSelection<-function(data,Type,variables,continuous, TD, NL,m=1,p=2,alpha_back=0.05,knots=-999){



  V<-length(variables) # number of total variables


  covB<-variables[order(1-continuous)]
  Ftd<-TD[order(1-continuous)]
  Fnl<-NL[order(1-continuous)]


  TDf<-rep(1,V)
  NLf<-c(rep(1,sum(continuous)),rep(0,sum(1-continuous)))
  nNLf<-sum(continuous)
  nNLfNOT=V-nNLf




  m1nl<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m1nl)<-0
  m2nl<-matrix(nrow=nNLf,ncol=nNLfNOT,0)
  m3nl<-cbind(m1nl,m2nl)
  m4nl<-matrix(nrow=V,ncol=V,c(rep(1,nNLf),rep(0,nNLfNOT)),byrow=T)
  mNL<-rbind(NLf,m3nl,m4nl)

  ## TD matrix
  m1td<-matrix(nrow=nNLf,ncol=V,1)
  m2td<-matrix(nrow=nNLfNOT,ncol=nNLf,1)
  m3td<-matrix(nrow=nNLfNOT,ncol=nNLfNOT,1)
  diag(m3td)<-0
  m4td<-cbind(m2td,m3td)
  m5td<-matrix(nrow=nNLf,ncol=nNLf,1)
  diag(m5td)<-0
  m6td<-matrix(nrow=nNLf,ncol=nNLfNOT,1)
  m7td<-cbind(m5td,m6td)
  mTD<-rbind(TDf,m1td,m4td,m7td)


  for (i in 1:V){
    if(Ftd[i]==1) mTD[,i]<-1
    else if(Ftd[i]==-1) mTD[,i]<-0

    if(Fnl[i]==1) mNL[,i]<-1
    else if(Fnl[i]==-1) mNL[,i]<-0
  }

  PLback<-rep(0,V+nNLf+1)
  DFback<-rep(-999,V+nNLf+1)

  for (i in 1:(V+nNLf+1)){
    res<-CoxFlex(data,Type,variables=covB,TD=mTD[i,],NL=mNL[i,],m,p,knots=-999)
    PLback[i]<-res$Partial_Log_Likelihood
    DFback[i]<-res$Number_of_parameters
  }


  Mind<-matrix(nrow=4,ncol=V,c(rep(1,2*V),rep(1,nNLf),rep(0,nNLfNOT),rep(2,nNLf),rep(1,nNLfNOT)),byrow=TRUE)
  for (i in 1:V){
    if(Ftd[i]==-1) Mind[2,i]<-0 # 2nd index TD
    if(Fnl[i]==-1) Mind[3,i]<-0 # 3rd index NL
  }
  Mind[4,]<-Mind[2,]+Mind[3,]


  pB<-rep(0,(V+nNLf)) # stores the test results of likelihood ratio tests for each model

  for (j in 1:(V+nNLf)){
    if(DFback[1]-DFback[j+1]!=0)
      pB[j]<-1-stats::pchisq(-2*(PLback[j+1]-PLback[1]),DFback[1]-DFback[j+1])
  }


  MpB<-matrix(nrow=5,ncol=V,0)

  MpB[4,1:nNLf]<-pB[1:nNLf]
  if(sum(1-continuous)!=0){
    MpB[2,(nNLf+1):V]<-pB[(nNLf+1):V]
  }
  MpB[5,1:nNLf]<-pB[(V+1):(V+nNLf)]


  a<-which.max(MpB)
  nc<-floor(a/5)+((a-5*floor(a/5))!=0)
  nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)

  indexrow<-1


  MatResPvalue<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1)
  MatResVariables<-matrix(nrow=2*length(variables)+sum(continuous),ncol=1)

  if(MpB[nr,nc]<alpha_back){
    mbase<-CoxFlex(data,Type,variables=covB,TD=mTD[1,],NL=mNL[1,],m=1,p=2,knots=-999)
  }

  while(MpB[nr,nc]>=alpha_back) {


    MatResPvalue[indexrow,1]<-MpB[nr,nc]
    MatResVariables[indexrow,1]<-paste(nr,nc,sep="")
    indexrow<-indexrow+1

    # print("indexrow",\n)

    if (nr==1) Mind[1,nc]<-0
    if (nr==2 | nr==5) Mind[2,nc]<-0
    if (nr==3 | nr==4) Mind[3,nc]<-0
    Mind[4,]<-Mind[2,]+Mind[3,]

    mbase<-CoxFlex(data,Type,variables=covB[Mind[1,]==1],TD=Mind[2,][Mind[1,]==1],NL=Mind[3,][Mind[1,]==1],m=1,p=2,knots=-999)

    MpB<-matrix(nrow=5,ncol=V,rep(0,V*5))


    for (k in 1:V){

      if (Mind[2,k]==1 & Mind[3,k]==0) {
        MindNew<-Mind
        if(Ftd[k]!=1) {
          MindNew[2,k]<-0
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[2,k]<-1-stats::pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)

        }
        else MpB[2,k]<-0
      }

      if (Mind[2,k]==1 & Mind[3,k]==1) {
        MindNew<-Mind
        if(Ftd[k]!=1){
          MindNew[2,k]<-0
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[5,k]<-1-stats::pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)

        }
        else MpB[5,k]<-0
      }

      if (Mind[3,k]==1 & Mind[2,k]==0) { #if NL + non-TD
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[3,k]<-1-stats::pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)

        }
        else MpB[3,k]<-0
      }

      if (Mind[3,k]==1 & Mind[2,k]==1) {
        MindNew<-Mind
        if(Fnl[k]!=1){
          MindNew[3,k]<-0
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[4,k]<-1-stats::pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)

        }
        else MpB[4,k]<-0
      }

      if (Mind[3,k]==0 & Mind[2,k]==0 & Mind[1,k]==1) {
        MindNew<-Mind
        if(Fnl[k]!=1 & Fnl[k]!=-1 & Ftd[k]!=1 & Ftd[k]!=-1){
          MindNew[1,k]<-0 # remove this variable
          mnew<-CoxFlex(data,Type,variables=covB[MindNew[1,]==1],TD=MindNew[2,][MindNew[1,]==1],NL=MindNew[3,][MindNew[1,]==1],m=1,p=2,knots=-999)
          MpB[1,k]<-1-stats::pchisq(-2*(mnew$Partial_Log_Likelihood-mbase$Partial_Log_Likelihood),mbase$Number_of_parameters-mnew$Number_of_parameters)

        }
        else MpB[1,k]<-0
      }
    }

    a<-which.max(MpB)
    nc<-floor(a/5)+((a-5*floor(a/5))!=0)
    nr<-a-5*floor(a/5)+5*((a-5*floor(a/5))==0)
  }

  list(final_model=mbase)
}
