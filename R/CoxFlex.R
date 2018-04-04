#' Estimate the TD and/or NL effect(s) of continuous covariate(s) in time-to-event analysis
#'
#' @param data A data frame in the long (interval) format with one line per unit of time.
#' @param Type A vector consisting the name of variables representing the start and stop of each
#'             time interval and event indicator. e.g. c("start","stop","event").
#' @param variables A vector consisting the name of varialbes that will be adjusted in the model.
#' @param TD A vector of binary indicators representing whether the corresponding variable has a time-dependent
#'           effect or not (1=Yes, 0=No)
#' @param NL A vector of binary indicators representing whether the corresponding variable has a non-linear
#'           effect or not (1=Yes, 0=No)
#' @param m The number of interior knots used in the regression B-spline
#' @param p The degree of the regression B-spline
#' @param knots Default value is -999
#' @return    Returns a list of the following items:
#' @return \item{PLL}{  partial loglikelihood of the final model}
#'         \item{NP}{  number of parameters estimated in the model}
#'         \item{NE}{  number of events in the dataset}
#'         \item{knots_covariates}{knots for splines of NL effects for corresponding covariates}
#'         \item{knots_time}{knots for splines of the TD effect}
#'         \item{coefficients_splines_NL}{estimated splines coefficients for NL effects}
#'         \item{coefficients_splines_TD}{estimated splines coefficients for TD effects}
#'         \item{variables}{name of the covariates in the order they are adjusted in the model}
#'         \item{coef}{estimated coefficients for covariates not having TD and NL effects}
#'         \item{sd}{standard errors of the corresponding estimated coefficients}
#'         \item{pval}{p-values for corresponding estimated effects}
#'
#' @export


CoxFlex<-function (data, Type, variables, TD, NL, m, p,knots) {

  res_principal<-last_prog(data, Type, variables, TD, NL, m, p,knots)

  TDtestP<-rep(-999,length(variables))
  TDtestN<-rep(-999,length(variables))
  TDtest<-rep(-999,length(variables))
  NLtestP<-rep(-999,length(variables))
  NLtestN<-rep(-999,length(variables))
  NLtest<-rep(-999,length(variables))

  for (kt in 1:length(variables)){
    if (TD[kt]==1){
      TDnew<-TD
      TDnew[kt]<-0
      resSEC<-last_prog(data, Type, variables, TDnew, NL, m, p,knots)
      TDtestP[kt]<-resSEC$Partial_Log_Likelihood
      TDtestN[kt]<-resSEC$Number_of_parameters
      TDtest[kt]<-1-stats::pchisq(-2*(TDtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-TDtestN[kt])
    }
    if (NL[kt]==1){
      NLnew<-NL
      NLnew[kt]<-0
      resSEC<-last_prog(data, Type, variables, TD, NLnew, m, p,knots)
      NLtestP[kt]<-resSEC$Partial_Log_Likelihood
      NLtestN[kt]<-resSEC$Number_of_parameters
      NLtest[kt]<-1-stats::pchisq(-2*(NLtestP[kt]-res_principal$Partial_Log_Likelihood),res_principal$Number_of_parameters-NLtestN[kt])
    }
  }



  Rescox<-matrix(ncol=6,nrow=sum(2*(NL+TD==2)+1*(NL+TD!=2)) )
  colnames(Rescox)<-c("","coef","exp(coef)","se(coef)","z","p")
  rownames(Rescox)<-rep(c(""),sum(2*(NL+TD==2)+1*(NL+TD!=2)))


  indexrescox<-0
  for (yu in 1:length(variables)){
    indexrescox<-indexrescox+1
    if (TD[yu]==0 & NL[yu]==0){
      Rescox[indexrescox,1]<-variables[yu]
      Rescox[indexrescox,2]<-round(as.numeric(res_principal$coefficients[yu]),3)
      Rescox[indexrescox,3]<-round(as.numeric(exp(res_principal$coefficients[yu])),3)
      Rescox[indexrescox,4]<-round(as.numeric(res_principal$Standard_Error[yu]),3)
      Rescox[indexrescox,5]<-round(as.numeric(res_principal$coefficients[yu]/res_principal$Standard_Error[yu]),3)
      Rescox[indexrescox,6]<-round(as.numeric(1-stats::pchisq((res_principal$coefficients[yu]/res_principal$Standard_Error[yu])^2,df=1)) ,3)
    }
    if (TD[yu]==0 & NL[yu]==1){
      Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(NLtest[yu],3)
    }
    if (TD[yu]==1 & NL[yu]==0){
      Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(TDtest[yu],3)
    }
    if (TD[yu]==1 & NL[yu]==1){
      Rescox[indexrescox,1]<-paste("NL(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(NLtest[yu],3)
      indexrescox<-indexrescox+1
      Rescox[indexrescox,1]<-paste("TD(",variables[yu],")",sep="")
      Rescox[indexrescox,2]<-"---"
      Rescox[indexrescox,3]<-"splines"
      Rescox[indexrescox,4]<-"---"
      Rescox[indexrescox,5]<-"---"
      Rescox[indexrescox,6]<-round(TDtest[yu],3)
    }
  }



  list(PLL = res_principal$Partial_Log_Likelihood,
       NP = res_principal$Number_of_parameters,
       NE=res_principal$Number_events ,
      # Number_knots  = res_principal$Number_knots, Degree_of_splines = res_principal$Degree_of_splines,
       knots_covariates = res_principal$knots_covariates,
       knots_time = res_principal$knots_time,
      # coefficients = res_principal$coefficients,
      # Standard_Error=res_principal$Standard_Error,
       coefficients_splines_NL = res_principal$coefficients_splines_NL,
       coefficients_splines_TD = res_principal$coefficients_splines_TD,
       variables=res_principal$variables,
       coef=suppressWarnings(as.numeric(Rescox[,2])),
       sd=suppressWarnings((as.numeric(Rescox[,2]))),
       pvalue=suppressWarnings(as.numeric(Rescox[,6])))



}
