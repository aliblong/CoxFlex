last_prog<-function (data, Type, variables, TD, NL, m, p,knots=-999) {





  if (length(Type)==2){Type2<-Type
  data$StartV0<-rep(0,dim(data)[1])
  Type<-c("StartV0",Type2[1],Type2[2])
  }

  i1 <- sum((NL + TD) == 0)
  i2 <- sum(((NL == 1) & (TD == 0)))
  i3 <- sum(((NL == 0) & (TD == 1)))
  i4 <- sum((NL + TD) == 2)
  nonpara <- TD + NL


  V <- length(variables)

  variablesNEW<-match(variables,names(data))
  TypeNEW<-match(Type,names(data))


  listeprobaquantile <- seq(1, m)/(m + 1)
  knotsNEW <- matrix(nrow = V + 1, ncol = p + 1 + m + p + 1)


  if (is.matrix(knots)==FALSE){ ## if not user defined knots
    if (is.numeric(knots)==TRUE & knots==-999) {
      for (i in 1:V) {
        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),
                           stats::quantile(data[, variables[i]], probs = listeprobaquantile),
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
      }
      knotsNEW[V + 1, ] <- c(rep(0, p + 1),
                             stats::quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile),
                             seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
    }
  }

  if (is.matrix(knots)==TRUE){ ## if user defined interior knots

    if (dim(knots)[1]!=(length(variables)+1) | dim(knots)[2]!=m) stop("Error Message: variable knots should be a matrix of dimention (length(variables)+1)*m")

    for (i in 1:V) {
      if (is.na(knots[i,1])==TRUE){  knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),
                                                        stats::quantile(data[, variables[i]], probs = listeprobaquantile),
                                                        seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
      ## if user defined knots for certain variable is NA, then use default knots setup
      } else {

        knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),
                           knots[i,],
                           seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        ## else use user defined interior knots
      }
    }

    if (is.na(knots[(V+1),1])==TRUE){  knotsNEW[(V+1), ] <- c(rep(0, p + 1),
                                                              stats::quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile),
                                                              seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
    } else {
      ## else use default knots for time
      knotsNEW[(V+1), ] <- c(rep(0, p + 1),
                             knots[(V+1),],
                             seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))

    }
  }


  data <- as.matrix(data)
  listeT <- c(0, sort(unique(data[data[, TypeNEW[3]] == 1, TypeNEW[2]])))
  ncol <- dim(data)[2]


  X <- split(data, data[, 1])
  matX <- sapply(X, DvlpMatrix, listeT = listeT, ncol = ncol, TypeNEW=TypeNEW)
  QWR <- do.call(rbind, matX)

  nbNL <- sum(NL)
  if (nbNL != 0) {
    for (i in 1:nbNL) {
      QWR <- cbind(QWR,  splines::splineDesign(knotsNEW[seq(1,V, 1)[NL == 1][i],], x=QWR[, variablesNEW[NL == 1][i]],ord=p+1)[,-1])
    }
  }

  QWR <- cbind(QWR, splines::splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]],ord=p+1))

  tt <- paste("modX<-survival::coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,",
              TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep = "")




  Nn <- 0
  nbpara <- sum((NL == 0 & TD == 0))

  if (nbpara != 0) {
    for (k in 1:nbpara) {
      tt <- paste(tt, "QWR[,", variablesNEW[NL == 0 & TD ==
                                              0][k], "]+", sep = "")
      Nn <- Nn + 1
    }
  }



  nbonlyNL <- sum((NL == 1 & TD == 0))
  if (nbonlyNL != 0) {
    onlyNL <- match(variablesNEW[NL == 1 & TD == 0], variablesNEW[NL ==
                                                                    1])
    for (k in 1:nbonlyNL) {
      covp<-paste( "QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p)+1):(dim(data)[2] + onlyNL[k]*(m+p))  , "]", sep = "")
      tt<-paste(tt, paste(c(covp,""), collapse= "+") )
      Nn <- Nn +m+p
    }
  }



  nbonlyTD <- sum((NL == 0 & TD == 1))
  if (nbonlyTD != 0) {
    for (k in 1:nbonlyTD) {
      flag<-dim(QWR)[2]+1
      QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD ==1][k]] * QWR[, (dim(data)[2] +(m+p)*nbNL
                                                                          +1): (dim(data)[2] +(m+p)*(nbNL+1)+1)])
      covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
      tt<-paste(tt, paste(c(covp,""), collapse= "+"))
      Nn <- Nn + m+p+1
    }#
  }

  tt2 <- tt


  nbNLTD <- sum((NL == 1 & TD == 1))
  NOTonlyNL <- match(variablesNEW[NL == 1 & TD == 1], variablesNEW[NL ==
                                                                     1])
  if (nbNLTD != 0) {
    for (k in 1:nbNLTD) {
      covp<-paste("QWR[,",(dim(data)[2] + (NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2] + NOTonlyNL[k]*(m+p)) , "]", sep = "" )
      tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
    }
    vrais <- c()
    modX <- survival::coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) -
                                             1))), method = "efron")
    vrais <- c(vrais, modX$loglik[2])


    mod<-modX
    tt2<-tt
    VV<-matrix(ncol=nbNLTD*(m+p+1),nrow=dim(QWR)[1])

    for (k in 1:nbNLTD){

      VV[,((1+m+p)*(k-1)+1):((1+m+p)*k)]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                              +(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]%*%mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))])
      covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
      tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))

    }
    modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
    vrais<-c(vrais,modX$loglik[2])




    diff<-1

    while(diff>0.00001){

      mod<-modX
      tt2<-tt
      VV<-matrix(ncol=nbNLTD*(m+p),nrow=dim(QWR)[1])
      for (k in 1:nbNLTD){
        VV[,((m+p)*(k-1)+1):((m+p)*(k))]<-QWR[,(dim(data)[2]+(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                          +(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]%*%mod$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])
        covp<-paste( "VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }

      modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2])



      mod<-modX
      tt2<-tt
      VV<-matrix(ncol=nbNLTD*(m+p+1),nrow=dim(QWR)[1])


      for (k in 1:nbNLTD){
        VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                  +(NOTonlyNL[k]-1)*(m+p)+1):(dim(data)[2]+NOTonlyNL[k]*(m+p))]%*%mod$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
        covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }



      modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2])
      diff<-abs(vrais[length(vrais)]-vrais[length(vrais)-2])


    }





  }    else {
    modX <- survival::coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) -
                                             1))), method = "efron")
    vrais <- c(modX$loglik[2])
  }
  rm(QWR, X, matX)

  MAT <- matrix(ncol = V, nrow = 1 + m + p + m + p + 1)
  if (i1 != 0) {
    for (j in 1:i1) {
      MAT[1, j] <- modX$coef[j]
    }
  }
  if (i2 != 0) {
    for (j in 1:i2) {
      MAT[2:(m + p + 1), i1 + j] <- modX$coef[(i1 + (j -
                                                       1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m +
                                                                            p)]
    }
  }
  if (i3 != 0) {
    for (j in 1:i3) {
      MAT[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- modX$coef[(i1 +
                                                                        i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 +
                                                                                                                     i2 * (m + p) + (j - 1) * (m + p + 1) + m + p +
                                                                                                                     1)]
    }
  }
  if (i4 != 0) {
    for (j in 1:i4) {
      MAT[2:(m + p + 1), i1 + i2 + i3 + j] <- mod$coef[(i1 +
                                                          i2 * (m + p) + i3 * (m + p + 1) + (j - 1) * (m +
                                                                                                         p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) +
                                                                                                                    (j - 1) * (m + p) + m + p)]
      MAT[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 +
            j] <- modX$coef[(i1 + i2 * (m + p) + i3 * (m +
                                                         p + 1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                (m + p) + i3 * (m + p + 1) + (j - 1) * (m + p +
                                                                                                                                          1) + m + p + 1)]
    }
  }


  MATse <- matrix(ncol = V, nrow = 1 + m + p + m + p + 1)
  if (i1 != 0) {
    for (j in 1:i1) {
      MATse[1, j] <- sqrt(diag(modX$var)[j])
    }
  }
  if (i2 != 0) {
    for (j in 1:i2) {
      MATse[2:(m + p + 1), i1 + j] <- sqrt(diag(modX$var)[(i1 + (j -
                                                                   1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m +
                                                                                        p)])
    }
  }
  if (i3 != 0) {
    for (j in 1:i3) {
      MATse[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- sqrt(diag(modX$var)[(i1 +
                                                                                    i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 +
                                                                                                                                 i2 * (m + p) + (j - 1) * (m + p + 1) + m + p +
                                                                                                                                 1)])
    }
  }
  if (i4 != 0) {
    for (j in 1:i4) {
      MATse[2:(m + p + 1), i1 + i2 + i3 + j] <- sqrt(diag(mod$var)[(i1 +
                                                                      i2 * (m + p) + i3 * (m + p + 1) + (j - 1) * (m +
                                                                                                                     p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) +
                                                                                                                                (j - 1) * (m + p) + m + p)])
      MATse[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 +
              j] <- sqrt(diag(modX$var)[(i1 + i2 * (m + p) + i3 * (m +
                                                                     p + 1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                            (m + p) + i3 * (m + p + 1) + (j - 1) * (m + p +
                                                                                                                                                      1) + m + p + 1)])
    }
  }


  var_order <- c(variablesNEW[NL == 0 & TD == 0], variablesNEW[NL ==
                                                                 1 & TD == 0], variablesNEW[NL == 0 & TD == 1], variablesNEW[NL ==
                                                                                                                               1 & TD == 1])

  coefficients<-MAT[1,match(variablesNEW,var_order)]
  names(coefficients)<-variables

  se_coef<-MATse[1,match(variablesNEW,var_order)]

  coefficients_splines_NL<-as.matrix(MAT[2:(m+p+1),match(variablesNEW,var_order)])
  coefficients_splines_NL<-rbind(rep(0,V),coefficients_splines_NL)
  coefficients_splines_NL[1,(NL==0)]<-NA
  colnames(coefficients_splines_NL)<-variables

  coefficients_splines_TD<-as.matrix(MAT[(m+p+2):(2*(m+p+1)),match(variablesNEW,var_order)])
  colnames(coefficients_splines_TD)<-variables



  knots_covariates<-knotsNEW[1:V,]
  if (V>1) {rownames(knots_covariates)<-variables}
  if (V>1) {knots_covariates[(NL==0),]<-rep(NA,p + 1 + m + p + 1)} else { if (NL==0) {knots_covariates<-rep(NA,p + 1 + m + p + 1)} }
  knots_time<-knotsNEW[V+1,]

  nEvents<-sum(data[,TypeNEW[3]]==1)

  rm(data, modX)
  gc()

  list(Partial_Log_Likelihood = vrais[length(vrais)], Number_of_parameters = nbpara +
         nbonlyNL * (m + p) + nbonlyTD * (1 + m + p) + nbNLTD *
         (m + p + m + p + 1), Number_events=nEvents, Number_knots = m, Degree_of_splines = p,
       knots_covariates = knots_covariates,
       knots_time = knots_time,
       coefficients = coefficients, Standard_Error=se_coef,coefficients_splines_NL = coefficients_splines_NL,coefficients_splines_TD = coefficients_splines_TD,variables=variables)


}
