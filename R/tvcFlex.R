#' Estimate the TD and/or NL effect(s) of time-varying continuous covariate(s) in time-to-event analysis
#'
#' @param data A data frame in the long (interval) format with one line per unit of time.
#' @param Type A vector consisting the name of variables representing the start/stop of each
#'             time interval, event indicator, and names of variable representing the TEL of each sparsely
#'             mearused time-varying covariates.  e.g. c("start","stop","event","tel").
#' @param variables A vector consisting the name of varialbes that will be adjusted in the model.
#' @param TD A vector of binary indicators representing whether the corresponding variable has a time-dependent
#'           effect or not (1=Yes, 0=No).
#' @param NL A vector of binary indicators representing whether the corresponding variable has a non-linear
#'           effect or not (1=Yes, 0=No).
#' @param TEL A vector of binary indicators representing whether to adjust the TEL effect of the corresponding variable
#'            (1=Yes, 0=No).
#' @param m The number of interior knots used in the regression B-spline
#' @param p The degree of the regression B-spline
#' @param knots Default value is -999, set the knots at the default place.
#' @return    Returns a list of the following items:
#' @return \item{PLL}{  partial loglikelihood of the final model}
#'         \item{NP}{  number of parameters estimated in the model}
#'         \item{NE}{  number of events in the dataset}
#'         \item{knots_covariates}{knots for splines of NL effects for corresponding covariates}
#'         \item{knots_time}{knots for splines of the TD effect}
#'         \item{knots_TEL}{knots for splines of the TEL effect}
#'         \item{variables}{name of the covariates in the order they are adjusted in the model}
#'         \item{coefficients_splines_NL}{estimated splines coefficients for NL effects}
#'         \item{coefficients_splines_TD}{estimated splines coefficients for TD effects}
#'         \item{coefficients_splines_TEL}{estimated splines coefficients for TEL effects}
#'         \item{coef}{estimated coefficients for covariates not having TD and NL effects}
#'         \item{sd}{standard errors of the corresponding estimated coefficients}
#'         \item{pval}{p-values for corresponding estimated effects}
#'@note Note that the TD and NL arguments have different meanings than that in the BackSelection()
#' @export




tvcFlex<-function (data, Type, variables, TD, NL,TEL, m, p,knots) {

  if (sum(TEL)==0){

    CoxFlex(data,Type,variables,TD,NL,m,p,knots)


  } else{
    res_principal<-last_prog_TEL(data, Type, variables, TD, NL,TEL ,m, p,knots)


    TDtestP<-rep(-999,length(variables))
    TDtestN<-rep(-999,length(variables))
    TDtest<-rep(-999,length(variables))
    TDdev<-rep(-999,length(variables))
    NLtestP<-rep(-999,length(variables))
    NLtestN<-rep(-999,length(variables))
    NLtest<-rep(-999,length(variables))
    NLdev<-rep(-999,length(variables))
    TELtestP<-rep(-999,length(variables))
    TELtestN<-rep(-999,length(variables))
    TELtest<-rep(-999,length(variables))
    TELdev<-rep(-999,length(variables))

    for (kt in 1:length(variables)){
      if (TD[kt]==1){
        TDnew<-TD
        TDnew[kt]<-0
        resSEC<-last_prog_TEL(data, Type, variables, TDnew, NL,TEL, m, p,knots)
        TDtestP[kt]<-resSEC$Partial_Log_Likelihood
        TDtestN[kt]<-resSEC$Number_of_parameters
        TDtest[kt]<-1-stats::pchisq(-2*(TDtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-TDtestN[kt])
        TDdev[kt]<--2*(TDtestP[kt]-res_principal$Partial_Log_Likelihood)
      }
      if (NL[kt]==1){
        NLnew<-NL
        NLnew[kt]<-0
        resSEC<-last_prog_TEL(data, Type, variables, TD, NLnew,TEL, m, p,knots)
        NLtestP[kt]<-resSEC$Partial_Log_Likelihood
        NLtestN[kt]<-resSEC$Number_of_parameters
        NLtest[kt]<-1-stats::pchisq(-2*(NLtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-NLtestN[kt])
        NLdev[kt]<--2*(NLtestP[kt]-res_principal$Partial_Log_Likelihood)
      }
      if (TEL[kt]==1){
        TELnew<-TEL
        TELnew[kt]<-0
        resSEC<-last_prog_TEL(data, Type, variables, TD, NL,TELnew, m, p,knots)
        TELtestP[kt]<-resSEC$Partial_Log_Likelihood
        TELtestN[kt]<-resSEC$Number_of_parameters
        TELtest[kt]<-1-stats::pchisq(-2*(TELtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-TELtestN[kt])
        TELdev[kt]<--2*(TELtestP[kt]-res_principal$Partial_Log_Likelihood)
      }
    }

    Rescox<-matrix(ncol=7,nrow=sum(3*(NL+TD+TEL==3)+2*(NL+TD+TEL==2)+1*(NL+TD+TEL!=2 & NL+TD+TEL!=3 )) )
    colnames(Rescox)<-c("","coef","exp(coef)","se(coef)","z","p","dev")
    rownames(Rescox)<-rep(c(""),sum(3*(NL+TD+TEL==3)+2*(NL+TD+TEL==2)+1*(NL+TD+TEL!=2 & NL+TD+TEL!=3 )))


    indexrescox<-0
    for (yu in 1:length(variables)){
      indexrescox<-indexrescox+1
      if (TD[yu]==0 & NL[yu]==0 & TEL[yu]==0){
        Rescox[indexrescox,1]<-variables[yu]
        Rescox[indexrescox,2]<-round(as.numeric(res_principal$coefficients[yu]),3)
        Rescox[indexrescox,3]<-round(as.numeric(exp(res_principal$coefficients[yu])),3)
        Rescox[indexrescox,4]<-round(as.numeric(res_principal$Standard_Error[yu]),3)
        Rescox[indexrescox,5]<-round(as.numeric(res_principal$coefficients[yu]/res_principal$Standard_Error[yu]),3)
        Rescox[indexrescox,6]<-round(as.numeric(1-stats::pchisq((res_principal$coefficients[yu]/res_principal$Standard_Error[yu])^2,df=1)) ,3)
        Rescox[indexrescox,7]<-res_principal$Partial_Log_Likelihood
      }
      if (TD[yu]==0 & NL[yu]==1 & TEL[yu]==0){
        Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(NLtest[yu],3)
        Rescox[indexrescox,7]<-round(NLdev[yu],4)
      }
      if (TD[yu]==1 & NL[yu]==0 & TEL[yu]==0){
        Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TDtest[yu],3)
        Rescox[indexrescox,7]<-round(TDdev[yu],4)
      }
      if (TD[yu]==0 & NL[yu]==0 & TEL[yu]==1){
        Rescox[indexrescox,1]<-paste("TEL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TELtest[yu],3)
        Rescox[indexrescox,7]<-round(TELdev[yu],4)
      }
      if (TD[yu]==1 & NL[yu]==1 & TEL[yu]==1){
        Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(NLtest[yu],3)
        Rescox[indexrescox,7]<-round(NLdev[yu],4)
        indexrescox<-indexrescox+1
        Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TDtest[yu],3)
        Rescox[indexrescox,7]<-round(TDdev[yu],4)
        indexrescox<-indexrescox+1
        Rescox[indexrescox,1]<-paste("TEL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TELtest[yu],3)
        Rescox[indexrescox,7]<-round(TELdev[yu],4)
      }
      if (TD[yu]==1 & NL[yu]==1 & TEL[yu]==0){
        Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(NLtest[yu],3)
        Rescox[indexrescox,7]<-round(NLdev[yu],4)
        indexrescox<-indexrescox+1
        Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TDtest[yu],3)
        Rescox[indexrescox,7]<-round(TDdev[yu],4)
      }
      if (TD[yu]==0 & NL[yu]==1 & TEL[yu]==1){
        Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(NLtest[yu],3)
        Rescox[indexrescox,7]<-round(NLdev[yu],4)
        indexrescox<-indexrescox+1
        Rescox[indexrescox,1]<-paste("TEL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TELtest[yu],3)
        Rescox[indexrescox,7]<-round(TELdev[yu],4)
      }
      if (TD[yu]==1 & NL[yu]==0 & TEL[yu]==1){
        Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TDtest[yu],3)
        Rescox[indexrescox,7]<-round(TDdev[yu],4)
        indexrescox<-indexrescox+1
        Rescox[indexrescox,1]<-paste("TEL(",variables[yu],")",sep="")
        Rescox[indexrescox,2]<-"---"
        Rescox[indexrescox,3]<-"splines"
        Rescox[indexrescox,4]<-"---"
        Rescox[indexrescox,5]<-"---"
        Rescox[indexrescox,6]<-round(TELtest[yu],3)
        Rescox[indexrescox,7]<-round(TELdev[yu],4)
      }
    }

    list(PLL= res_principal$Partial_Log_Likelihood,
         NP= res_principal$Number_of_parameters,
         NE=res_principal$Number_events,
         knots_covariates = res_principal$knots_covariates,
         knots_time = res_principal$knots_time,
         knots_TEL = res_principal$knots_TEL,
         variables=res_principal$variables,
       #  coefficients = res_principal$coefficients, Standard_Error=res_principal$Standard_Error,
         coefficients_splines_NL = res_principal$coefficients_splines_NL,
         coefficients_splines_TD = res_principal$coefficients_splines_TD,
         coefficients_splines_TEL = res_principal$coefficients_splines_TEL,
         coef=suppressWarnings(as.numeric(Rescox[,2])),
         sd=suppressWarnings((as.numeric(Rescox[,4]))),
         pvalue=suppressWarnings(as.numeric(Rescox[,6])))
        # dev=suppressWarnings(as.numeric(Rescox[,7])))

  }

}
