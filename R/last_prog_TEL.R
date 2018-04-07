last_prog_TEL<-function (data, Type, variables, TD, NL,TEL, m, p,knots=-999) {


  if(sum(TEL)==0) { # if there is no TEL effect

    last_prog(data=data, Type,variables,TD,NL,m,p,knots = -999)


  } else{   ## if there are TEL effect

    i1 <- sum((NL+TD+TEL) == 0)    ##number of non-nl non-td non-TEL variables
    i2 <- sum(((NL == 1) & (TD == 0) & (TEL==0))) ##number of NL non-td non-TEL variables
    i3 <- sum(((NL == 0) & (TD == 1) & (TEL==0)))  ##number of non-NL td non-TEL variables
    i4 <- sum((NL + TD) == 2 & (TEL==0)) ###number of NL TD non-TEL variables
    i5 <- sum(((NL == 0) & (TD == 0) & (TEL==1))) ##number of non-NL non-td TEL variables
    i6 <- sum(((NL == 1) & (TD == 0) & (TEL==1))) ##number of NL non-td TEL variables
    i7 <- sum(((NL == 0) & (TD == 1) & (TEL==1))) ##number of non-NL td TEL variables
    i8 <- sum(((NL == 1) & (TD == 1) & (TEL==1))) ##number of NL td TEL variables

    V <- length(variables) # number of total covariates

    nTEL=sum(TEL) # number of variables having TEL effect

    variablesNEW<-match(variables,names(data))  # stores the postition of each variable in the newdata set
    TypeNEW<-match(Type,names(data)) # stored the position of the Type varialbes in the data set


    listeprobaquantile <- seq(1, m)/(m + 1) ## place the interior knots according to quantiles
    knotsNEW <- matrix(nrow = V + 1+nTEL, ncol = p + 1 + m + p + 1) # stores the total knots (2(p+1)+m) for all the varialbes in the data set
    # and the knots for the time in all (v+1)row, add knots for TEL on (V+nTEL) row

    if (is.matrix(knots)==FALSE){ ## if not user defined knots
      if (is.numeric(knots)==TRUE & knots==-999) { ## if the knots is set as default
        for (i in 1:V) {
          knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1), ## place the first p+1 exterior knots at the min(variables[i])
                             stats::quantile(data[, variables[i]], probs = listeprobaquantile), ## place the interior knots at quantile
                             seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1)) ### place the last p+1 exterior knots equally spaced between (max(variables[i]), ...+p)
        }
        knotsNEW[V + 1, ] <- c(rep(0, p + 1), ## place the exterior knots for time at 0
                               stats::quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile), ### place the interior knots at quantile of event times
                               seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
        # place exterior knots equally spaced between max(event time) and max(event time)+p

        for(i in 1:nTEL){
          knotsNEW[(V+1+i), ] <- c(rep(0, p + 1), ## place the exterior knots for time at 0
                                   stats::quantile(data[, TypeNEW[3+i]], probs = 0.25), ### place the interior knots at quantile of event times
                                   seq(max(data[, TypeNEW[3+i]]), max(data[ ,TypeNEW[3+i]]) + p, 1))
          # place exterior knots equally spaced between max(event time) and max(event time)+p
        }
      }
    }

    if (is.matrix(knots)==TRUE){ ## if user defined interior knots

      if (dim(knots)[1]!=(length(variables)+2) | dim(knots)[2]!=m) stop("Error Message: variable knots should be a matrix of dimention (length(variables)+2)*m")

      for (i in 1:V) {
        if (is.na(knots[i,1])==TRUE){  knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),
                                                          stats::quantile(data[, variables[i]], probs = listeprobaquantile),
                                                          seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
        ## if user defined knots for certain variable is NA, then use default knots setup
        } else {

          knotsNEW[i, ] <- c(rep(min(data[, variables[i]]), p + 1),
                             knots[i,],
                             seq(max(data[, variables[i]]), max(data[, variables[i]])+p, 1))
             }
      }

      if (is.na(knots[(V+1),1])==TRUE){  knotsNEW[(V+1), ] <- c(rep(0, p + 1),
                                                                stats::quantile(data[data[, TypeNEW[3]] == 1, TypeNEW[2]], probs = listeprobaquantile),
                                                                seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))
      # if user defined knots for time is NA, use default
      } else {
        ## else use default knots for time
        knotsNEW[(V+1), ] <- c(rep(0, p + 1),
                               knots[(V+1),],
                               seq(max(data[data[,TypeNEW[3]] == 1, TypeNEW[2]]), max(data[data[, TypeNEW[3]] ==1, TypeNEW[2]]) + p, 1))

      }
      if (is.na(knots[(V+2),1])==TRUE){  knotsNEW[(V+2), ] <- c(rep(0, p + 1),
                                                                stats::quantile(data[, TypeNEW[4]], probs = 0.25),
                                                                seq(max(data[, TypeNEW[4]]), max(data[, TypeNEW[4]]) + p, 1))
      # if user defined knots for TEL is NA, use default
      } else {
        ## else use default knots for time
        knotsNEW[(V+2), ] <- c(rep(0, p + 1),
                               knots[(V+2),],
                               seq(max(data[, TypeNEW[4]]), max(data[, TypeNEW[4]]) + p, 1))

      }

    }



    data <- as.matrix(data)
    listeT <- c(0, sort(unique(data[data[, TypeNEW[3]] == 1, TypeNEW[2]]))) # sort the unique event time and stored in listeT
    ncol <- dim(data)[2] # number of coloum


    X <- split(data, data[, 1]) ##lists store full data for each individual respectively
    matX <- sapply(X, DvlpMatrix, listeT = listeT, ncol = ncol,TypeNEW=TypeNEW) ## transform the data structure
    QWR <- do.call(rbind, matX) #combine the lists to a full data matrix

    nbNL <- sum(NL) # number of NL effects
    if (nbNL != 0) {  # generate the spline values for NL effect, since the constraint a=0 was added to all the first spline
      for (i in 1:nbNL) {
        QWR <- cbind(QWR,  splines::splineDesign(knotsNEW[seq(1,V, 1)[NL == 1][i],], x=QWR[, variablesNEW[NL == 1][i]],ord=p+1)[,-1])
      }
    }


    # spline value for time
    QWR <- cbind(QWR, splines::splineDesign(knotsNEW[V+1,], x=QWR[, TypeNEW[2]],ord=p+1))


    for(j in 1:nTEL){
      # spline value for TEL
      QWR <- cbind(QWR,splines::splineDesign(knotsNEW[(V+1+j),], x=QWR[, TypeNEW[3+j]],ord=p+1))
    }


    tt <- paste("modX<-survival::coxph(Surv(QWR[,", TypeNEW[1], "],QWR[,",
                TypeNEW[2], "],QWR[,", TypeNEW[3], "])~", sep = "")

    ############## if there is none or only 1 special effect #####
    Nn <- 0
    nbpara <- sum((NL == 0 & TD == 0 & TEL==0 )) # number of parameters have neither NL nor TD effects nor TEL effects
    if (nbpara != 0) { # if there are variables have neither NL nor TD nor TEL effects
      for (k in 1:nbpara) {
        tt <- paste(tt, "QWR[,", variablesNEW[NL == 0 & TD ==
                                                0 & TEL==0][k], "]+", sep = "")
        Nn <- Nn + 1
      }
    }##Nn=number NL=TD=TEL=0


    nbonlyNL <- sum((NL == 1 & TD == 0 & TEL==0)) # only NL no TD no TEL
    if (nbonlyNL != 0) {
      onlyNL <- match(variablesNEW[NL == 1 & TD == 0 & TEL==0], variablesNEW[NL ==
                                                                               1]) # index the order of the varialbes only NL no TD no TEL effects among those
      # have NL effects, i.e., the order of the spline coloum
      for (k in 1:nbonlyNL) {
        covp<-paste( "QWR[,", (dim(data)[2] + (onlyNL[k]-1)*(m+p)+1):(dim(data)[2] + onlyNL[k]*(m+p))  , "]", sep = "")
        tt<-paste(tt, paste(c(covp,""), collapse= "+") )
        Nn <- Nn +m+p
      }
    }##Nn=number NL=TD=TEL=0 +(TD=TEL=0, NL=1)
    nbonlyTD <- sum((NL == 0 & TD == 1 & TEL==0)) # no NL no TEL only TD
    if (nbonlyTD != 0) {
      for (k in 1:nbonlyTD) {
        flag<-dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD ==1& TEL==0][k]] * QWR[, (dim(data)[2] +(m+p)*nbNL
                                                                                    +1): (dim(data)[2] +(m+p)*(nbNL+1)+1)])
        covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
        tt<-paste(tt, paste(c(covp,""), collapse= "+"))  #### only the TD spline base combined with covariate values
        # will enter the model
        Nn <- Nn + m+p+1

      }##Nn=number NL=TD=TEL=0 +(TD=TEL=0, NL=1)*nbonlyNL +(NL=TEL=0,TD=1)*nbonlyTD
    }
    nbonlyTEL <- sum((NL == 0 & TD == 0 & TEL==1)) # no NL no TD only TEL
    if (nbonlyTEL != 0) {
      onlyTEL <- match(variablesNEW[NL == 0 & TD == 0 & TEL==1], variablesNEW[TEL ==
                                                                                1]) # index the order of the varialbes only TEL no TD no NL effects among those

      for (k in 1:nbonlyTEL) {
        flag<-dim(QWR)[2]+1
        QWR <- cbind(QWR, QWR[, variablesNEW[NL == 0 & TD == 0 & TEL==1][k]] * QWR[, (dim(data)[2] +(m+p)*(nbNL+1)+1
                                                                                      +(onlyTEL[k]-1)*(m+p+1)+1):(dim(data)[2] +(m+p)*(nbNL+1)+1+onlyTEL[k]*(1+m+p))])
        covp<-paste("QWR[,", flag: dim(QWR)[2], "]", sep = "" )
        tt<-paste(tt, paste(c(covp,""), collapse= "+"))  # add TEL spline base multiplied by the covariates value
        Nn <- Nn +1+m+p

      }##Nn=number NL=TD=TEL=0 +(TD=TEL=0, NL=1)*nbonlyNL +(NL=TEL=0,TD=1)*nbonlyTD +(NL=TD=0,TEL=1)*nbonlyTEL
    }
    tt2 <- tt ## tt mark the non or only 1 effect estimates

    ### having more than one special effects###

    nbNLTDTEL<-sum((NL == 1 & TD == 1 & TEL==1)) # both NL and TD and TEL
    pNLTDTEL <- match(variablesNEW[NL == 1 & TD == 1 & TEL==1], variablesNEW[NL == 1]) # mark these variable

    nbNLTD <- sum((NL == 1 & TD == 1 & TEL==0)) # both NL and TD no TEL
    pNLTD <- match(variablesNEW[NL == 1 & TD == 1 & TEL==0], variablesNEW[NL ==1])  # mark these variable

    nbNLTEL <- sum((NL == 1 & TD == 0 & TEL==1)) # both NL and TEL no TD
    pNLTEL <- match(variablesNEW[NL == 1 & TD == 0 & TEL==1], variablesNEW[NL ==1])  # mark these variable

    nbTDTEL <- sum((NL == 0 & TD == 1 & TEL==1)) # both TD and TEL no NL
    pTDTEL <- match(variablesNEW[NL == 0 & TD == 1 & TEL==1], variablesNEW[TD ==1])  # mark these variable


    posTELNLTD <- match(variablesNEW[NL == 1 & TD == 1 & TEL==1], variablesNEW[TEL ==1])  # mark these variable

    posTELTD <- match(variablesNEW[NL == 0 & TD == 1 & TEL==1], variablesNEW[TEL ==1])  # mark these variable

    posTELNL <- match(variablesNEW[NL == 1 & TD == 0 & TEL==1], variablesNEW[TEL ==1])  # mark these variable


    vrais <- c()

    if(nbNLTDTEL !=0){ # having 3 effects

      ## starting with estimating NL effect first#####

      for(k in 1:nbNLTDTEL){ # for the covariates having 3 effects, estimate NL first
        covp<-paste("QWR[,", (dim(data)[2] + (pNLTDTEL[k] -
                                                1) * (m + p) +1) : (dim(data)[2]+ pNLTDTEL[k] * (m + p)), "]", sep = "" )
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }

      if (nbNLTD != 0) { # if there is both NL and TD but no TEL effect
        for (k in 1:nbNLTD) {
          covp<-paste("QWR[,", (dim(data)[2] + (pNLTD[k] -
                                                  1) * (m + p) +1) : (dim(data)[2]+ (pNLTD[k])* (m + p)), "]", sep = "" )
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }
      }


      if (nbNLTEL != 0) { # if there is both NL and TEL but no TD effect
        for (k in 1:nbNLTEL) {

          covp<-paste("QWR[,", (dim(data)[2] + (pNLTEL[k] -
                                                  1) * (m + p) +1) : (dim(data)[2]+ (pNLTEL[k]) * (m + p)), "]", sep = "" )
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }
      }


      modX <- survival::coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) -
                                               1))), method = "efron")
      vrais <- c(vrais, modX$loglik[2]) # estimate all the NL effects first

      modNL<-modX
      tt2<-tt
      VV<-matrix(ncol=(nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])

      ##Estimate TD effect conditional on NL effect

      # estimate TD effects (TDNLTEL+TDNL+TDTEL)
      for (k in 1:nbNLTDTEL){

        kal<-QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]

        VV[,((m+p+1)*(k-1)+1): (k*(m+p+1))]<-as.vector(kal)*QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+1+m+p)]
        covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }
      if (nbNLTD !=0){
        for (k in 1:nbNLTD){
          VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                        +(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))])
          covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }#spline for TD*NL spline of the variables having NLTD*previouly estimated corresponding coef
      }

      if (nbTDTEL !=0){
        for (k in 1:nbTDTEL){## only TD spline and the variables having these effect, no TEL yet
          VV[,((1+m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTD+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]

          covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
        }#spline for TD*the covariates having TDTEL effects
      }

      modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2]) # estimate the TD effects conditional on the NL effects

      ############################################
      # estimate TEL effects conditional on NL TD(TDNLTEL+NLTEL+TDTEL)
      modTD<-modX
      tt2<-tt ## VV stores the TEL spline
      VV<-matrix(ncol=(nbNLTDTEL+nbNLTEL+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
      for (k in 1:nbNLTDTEL){
        kal<-QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]
        kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))]%*%modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]
        VV[,((m+p+1)*(k-1)+1): (k*(m+p+1))]<-as.vector(kal)*as.vector(kalTD)*QWR[,(dim(data)[2]+(m+p)*nbNL+1+m+p+(posTELNLTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+1+m+p+posTELNLTD[k]*(m+p+1))]
        covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep = "")
        tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
      }
      if (nbNLTEL !=0){
        for (k in 1:nbNLTEL){
          VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1)+posTELNL[k]*(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                                                          +(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))])
          covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
        }
      }


      if (nbTDTEL !=0){
        for (k in 1:nbTDTEL){

          VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTEL+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1)+posTELTD[k]*(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                                                                          +nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1))]%*%modTD$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]])

          covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTEL+k-1)+1):((k+nbNLTEL+nbNLTDTEL)*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
        }
      }


      modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
      vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects


      diff<-1

      while(diff>0.0001){ ## ACE algorithm until converge

        ##########################################################################
        # estimate NL effect first conditional on previouly estimated TD, TEL effect
        modTEL<-modX
        tt2<-tt
        VV<-matrix(ncol=(nbNLTDTEL+nbNLTD+nbNLTEL)*(m+p),nrow=dim(QWR)[1])
        for (k in 1:nbNLTDTEL){
          kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)]%*%modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]
          kalTEL<-QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1)+(posTELNLTD[k])*(m+p+1))]%*%modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]
          VV[,((m+p)*(k-1)+1):(k*(m+p))]<-as.vector(kalTD)*as.vector(kalTEL)*QWR[,(dim(data)[2] + (pNLTDTEL[k]-1) * (m + p) +1) : (dim(data)[2]+ pNLTDTEL[k] * (m + p))]
          covp<-paste( "VV[,", ((m+p)*(k-1)+1):(k*(m+p)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }
        if(nbNLTD !=0){
          for (k in 1:nbNLTD){
            VV[,((m+p)*(nbNLTDTEL+k-1)+1):((m+p)*(nbNLTDTEL+k))]<-QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                          +(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]%*%modTD$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))])
            covp<-paste( "VV[,", ((m+p)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))

          }#NL spline base*TD base *TD coef
        }

        if(nbNLTEL !=0){

          for (k in 1:nbNLTEL){
            VV[,((m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((m+p)*(nbNLTDTEL+nbNLTD+k))]<-QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                          +nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNL[k]*(m+p+1))]%*%modTEL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+k)*(m+p+1))])

            covp<-paste( "VV[,", ((m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )

          }#NL spline base*TEL base *TEL coef

        }

        modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2])

        ########################################################
        # estimate TD effect conditional on the NL and TEL effects

        modNL<-modX
        tt2<-tt
        VV<-matrix(ncol=(nbNLTDTEL+nbNLTD+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
        for (k in 1:nbNLTDTEL){

          kal<-QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]
          kalTEL<-QWR[,(dim(data)[2]+nbNL*(m+p)+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(m+p+1)+posTELNLTD[k]*(m+p+1))]%*%modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]
          VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))]<-as.vector(kal)*as.vector(kalTEL)*QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+1+m+p)]
          covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }


        if(nbNLTD !=0){
          for (k in 1:nbNLTD){
            VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                          +(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTDTEL+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+k)*(m+p))])
            covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))

          }
        }
        if(nbTDTEL !=0){
          for (k in 1:nbTDTEL){
            VV[,((1+m+p)*(nbNLTDTEL+nbNLTD+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTD+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                        +nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(posTELTD[k]+1)*(m+p+1))]%*%modTEL$coef[(Nn+(nbNLTDTEL+nbNLTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+nbNLTEL+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]])

            covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTD+k-1)+1):((k+nbNLTD+nbNLTDTEL)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )

          }
        }

        modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2])

        ###############################################################################
        # estimate TEL effects conditional on the NL and TD effects(TDNLTEL+NLTEL+TDTEL)


        modTD<-modX
        tt2<-tt
        VV<-matrix(ncol=(nbNLTDTEL+nbNLTEL+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
        for (k in 1:nbNLTDTEL){
          kal<-QWR[,(dim(data)[2]+(pNLTDTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTDTEL[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+k*(m+p))]
          # NL base multiply by the NL coefficients
          kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)]%*%modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+k*(m+p+1))]
          VV[,((m+p+1)*(k-1)+1):(k*(m+p+1))]<-as.vector(kal)*as.vector(kalTD)*QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNLTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELNLTD[k]+1)*(m+p+1))] #TEL splinie base multiply kal
          covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):(k*(m+p+1)), "]", sep = "")
          tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
        }
        if (nbNLTEL !=0){
          for (k in 1:nbNLTEL){
            VV[,((1+m+p)*(nbNLTDTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                                                        +(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p))])
            covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+k-1)+1):((k+nbNLTDTEL)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
        }
        if (nbTDTEL !=0){
          for (k in 1:nbTDTEL){
            kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)]%*%modTD$coef[(Nn+(nbNLTDTEL+nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTDTEL+nbNLTD+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]

            VV[,((1+m+p)*(nbNLTDTEL+nbNLTEL+k-1)+1):((1+m+p)*(nbNLTDTEL+nbNLTEL+k))]<-as.vector(kalTD)*QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]

            covp<-paste( "VV[,", ((m+p+1)*(nbNLTDTEL+nbNLTEL+k-1)+1):((k+nbNLTEL+nbNLTDTEL)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }
        }


        modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects


        diff<-vrais[length(vrais)]-vrais[length(vrais)-3]
        print(diff)
      }


    } else if(nbNLTD+nbNLTEL+nbTDTEL !=0){ # having only two of the three effects
      #####################################
      ## estimate only the NL effect first if there is any#


      if(nbNLTD+nbNLTEL !=0){ #if there is a NL effect

        if (nbNLTD != 0) { # if there is both NL and TD but no TEL effect
          for (k in 1:nbNLTD) {
            covp<-paste("QWR[,", (dim(data)[2] + (pNLTD[k] -
                                                    1) * (m + p) +1) : (dim(data)[2]+ (pNLTD[k])* (m + p)), "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }

        if(nbNLTEL !=0){ # if there is NL TEL, no TD
          for (k in 1:nbNLTEL) {
            covp<-paste("QWR[,", (dim(data)[2] + (pNLTEL[k] -
                                                    1) * (m + p) +1) : (dim(data)[2]+ (pNLTEL[k]) * (m + p)), "]", sep = "" )
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }
        }

        modX <- survival::coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) -
                                                 1))), method = "efron")
        vrais <- c(vrais, modX$loglik[2])


        modNL<-modX
        tt2<-tt
      }
      ##########################################
      ## estimate TD effect conditional on the NL
      if(nbNLTD+nbTDTEL !=0){ #if there is a TD effect
        VV<-matrix(ncol=(nbNLTD+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
        if(nbNLTD!=0){
          for (k in 1:nbNLTD){
            VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                      +(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
            covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
          }# TD spline* NL spline*coef NL
        }

        if(nbTDTEL!=0){
          for (k in 1:nbTDTEL){
            VV[,((1+m+p)*(nbNLTD+k-1)+1):((1+m+p)*(nbNLTD+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]

            covp<-paste( "VV[,", ((m+p+1)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }#TD spline*corresponding variable
        }
        modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2])

        modTD<-modX
        tt2<-tt
      }

      ############################
      #estimate TEL conditional on TD NL

      if(nbNLTEL+nbTDTEL !=0){ #if there is a TEL effect
        VV<-matrix(ncol=(nbNLTEL+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])

        if (nbNLTEL !=0){
          for (k in 1:nbNLTEL){
            VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                                    +(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))])
            covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }#TEL spline*(sum(NL spline*corresponding covariates))
        }


        if (nbTDTEL !=0){
          for (k in 1:nbTDTEL){
            kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)]%*%modTD$coef[(Nn+(nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTD+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]
            VV[,((1+m+p)*(nbNLTEL+k-1)+1):((1+m+p)*(nbNLTEL+k))]<-as.vector(kalTD)*QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]

            covp<-paste( "VV[,", ((m+p+1)*(nbNLTEL+k-1)+1):((k+nbNLTEL)*(m+p+1)), "]", sep = "")
            tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
          }#TEL spline base*corresponding covariates*sum(TD spline base*corresponding coef)
        }

        modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
        vrais<-c(vrais,modX$loglik[2]) # estimate the TEL effects conditional on the TD and NL effects

        modTEL<-modX
        tt2<-tt
      }

      diff<-1
      #w<-1
      #w<-0
      while(diff>0.0001){ ## ACE algorithm until converge
        #############################
        #estimate NL effect conditional on TD, TEL effects


        if(nbNLTD+nbNLTEL !=0){ #if there is a NL effect

          VV<-matrix(ncol=(nbNLTD+nbNLTEL)*(m+p),nrow=dim(QWR)[1])

          if(nbNLTD!=0){
            for (k in 1:nbNLTD){
              VV[,((m+p)*(k-1)+1):((m+p)*(k))]<-QWR[,(dim(data)[2]+(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                        +(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]%*%modTD$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])
              covp<-paste( "VV[,", ((m+p)*(k-1)+1):((k)*(m+p)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))

            }#NL spline*(sum(TD spline*TD coef)
          }

          if(nbNLTEL!=0){
            for (k in 1:nbNLTEL){
              VV[,((m+p)*(nbNLTD+k-1)+1):((m+p)*(nbNLTD+k))]<-QWR[,(dim(data)[2]+(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                        +nbNL*(m+p)+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(posTELNL[k]+1)*(m+p+1))]%*%modTEL$coef[(Nn+(k-1)*(m+p+1)+1):(Nn+(k)*(m+p+1))])

              covp<-paste( "VV[,", ((m+p)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            }#NL spline*(sum(TEL spline*TEL coef))
          }
          modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])

          modNL<-modX
          tt2<-tt
        }
        ##########################################################
        # estimate TD conditional on NL,TEL

        if(nbNLTD+nbTDTEL !=0){ #if there is a TD effect
          VV<-matrix(ncol=(nbNLTD+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
          if(nbNLTD!=0){
            for (k in 1:nbNLTD){
              VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                        +(pNLTD[k]-1)*(m+p)+1):(dim(data)[2]+pNLTD[k]*(m+p))]%*%modNL$coef[(Nn+(k-1)*(m+p)+1):(Nn+(k)*(m+p))])
              covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+"))
            }#TD splines*(sum(NL splines*coef NL))
          }

          if(nbTDTEL!=0){
            for (k in 1:nbTDTEL){
              VV[,((1+m+p)*(nbNLTD+k-1)+1):((1+m+p)*(nbNLTD+k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+1):(dim(data)[2]+(m+p)*nbNL+(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                      +nbNL*(m+p)+m+p+1+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+nbNL*(m+p)+(posTELTD[k]+1)*(m+p+1))]%*%modTEL$coef[(Nn+(nbNLTEL+k-1)*(m+p+1)+1):(Nn+(nbNLTEL+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]])

              covp<-paste( "VV[,", ((m+p+1)*(nbNLTD+k-1)+1):((k+nbNLTD)*(m+p+1)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            }#TD spline*sum(TEL spline*TEL coef)*corresponding variable
          }

          modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])

          modTD<-modX
          tt2<-tt

        }

        ###############################################################################
        # estimate TEL effects conditional on the NL and TD effects (NLTEL+TDTEL)

        if(nbNLTEL+nbTDTEL !=0){ #if there is a TEL effect
          VV<-matrix(ncol=(nbNLTEL+nbTDTEL)*(m+p+1),nrow=dim(QWR)[1])
          if (nbNLTEL !=0){
            for (k in 1:nbNLTEL){
              VV[,((1+m+p)*(k-1)+1):((1+m+p)*(k))]<-QWR[,(dim(data)[2]+(m+p)*nbNL+m+p+1+(posTELNL[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELNL[k]+1)*(m+p+1))]*as.vector(QWR[,(dim(data)[2]
                                                                                                                                                                                      +(pNLTEL[k]-1)*(m+p)+1):(dim(data)[2]+pNLTEL[k]*(m+p))]%*%modNL$coef[(Nn+(nbNLTD+k-1)*(m+p)+1):(Nn+(nbNLTD+k)*(m+p))])
              covp<-paste( "VV[,", ((m+p+1)*(k-1)+1):((k)*(m+p+1)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )
            }# TEL spline*sum(NL spline*NL coef)
          }
          if (nbTDTEL !=0){
            for (k in 1:nbTDTEL){

              kalTD<-QWR[,(dim(data)[2]+nbNL*(m+p)+1):(dim(data)[2]+nbNL*(m+p)+m+p+1)]%*%modTD$coef[(Nn+(nbNLTD+k-1)*(m+p+1)+1):(Nn+(nbNLTD+k)*(m+p+1))]*QWR[,variablesNEW[NL == 0 & TD == 1 & TEL==1][k]]

              VV[,((1+m+p)*(nbNLTEL+k-1)+1):((1+m+p)*(nbNLTEL+k))]<-as.vector(kalTD)*QWR[,(dim(data)[2]+(m+p)*nbNL+(m+p+1)+(posTELTD[k]-1)*(m+p+1)+1):(dim(data)[2]+(m+p)*nbNL+(posTELTD[k]+1)*(m+p+1))]

              covp<-paste( "VV[,", ((m+p+1)*(nbNLTEL+k-1)+1):((k+nbNLTEL)*(m+p+1)), "]", sep = "")
              tt2<-paste(tt2, paste(c(covp,""), collapse= "+") )

            }#TEL spline*sum(TD spline* TD coef)* corresponding covariates
          }


          modX<-survival::coxph(eval(parse(text=substr(tt2,13,nchar(tt2)-1))),method="efron")
          vrais<-c(vrais,modX$loglik[2])
          modTEL<-modX
          tt2<-tt
        }

        if((i4!=0&i6!=0)|(i4!=0&i7!=0)|(i6!=0&i7!=0)){
          diff<-vrais[length(vrais)]-vrais[length(vrais)-3]}
        else{
          diff<-vrais[length(vrais)]-vrais[length(vrais)-2]
        }
        ### since there will only be two of the effects, so compare every two models
        # diff<-vrais[length(vrais)]-vrais[length(vrais)-2]
        print(diff)

      }

    }    else {
      modX <- survival::coxph(eval(parse(text = substr(tt2, 13, nchar(tt2) -
                                               1))), method = "efron")
      vrais <- c(modX$loglik[2])
    }
    rm(QWR,X,matX)

    MAT <- matrix(ncol = V, nrow = 1 + m + p + m + p+1+m+p+1) # store the estimated parameter
    if (i1 != 0) { #if no NL and TD no TEL
      for (j in 1:i1) { # first row store estimated coefficients for those having no NL TD TEL effects
        MAT[1, j] <- modX$coef[j]
      }
    }
    if (i2 != 0) {# only NL no TD no TEL
      for (j in 1:i2) { # 2nd to m+p+1 row stores the m+p coefficients for thoes only having NL
        # effects, each column represents one variable
        MAT[2:(m + p + 1), i1 + j] <- modX$coef[(i1 + (j -
                                                         1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m +
                                                                              p)]
      }
    }
    if (i3 != 0) { # only TD no NL no TEL
      for (j in 1:i3) { # row m+p+2 to 2*(m+p+1) stores TD effects
        MAT[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- modX$coef[(i1 +
                                                                          i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 +
                                                                                                                       i2 * (m + p) + (j - 1) * (m + p + 1) + m + p +1)]

      }
    }

    if (i5 != 0) { # only TEL no NL no TD
      for (j in 1:i5) { # row m+p+2 to 2*(m+p+1) stores TD effects
        MAT[(2*m + 2*p + 3):(3 * m + 3 * p + 3), i1 + i2 + i3+j] <- modX$coef[(i1 +
                                                                                 i2 * (m + p) + i3 * (m + p + 1) +(j - 1) * (m + p + 1)+ 1):(i1 +  i2 * (m + p) +i3 * (m + p + 1)+ (j - 1) * (m + p + 1) + m + p + 1)]
      }
    }

    if(i8!=0){ #NL+TD+TEL
      for (j in 1:i8) {
        MAT[2:(m + p + 1), i1 + i2 + i3+i5 + j] <- modNL$coef[(i1 +
                                                                 i2 * (m + p) + i3 * (m + p + 1)+i5*(m+p+1) + (j - 1) * (m + p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) + i5*(m+p+1)+(j - 1) * (m + p) + m + p)]
        MAT[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5+
              j] <- modTD$coef[(i1 + i2 * (m + p) + i3 * (m +  p + 1)+i5*(m+p+1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                 (m + p) + i3 * (m + p + 1) + i5*(m+p+1) + (j - 1) * (m + p + 1) + m + p + 1)]
        MAT[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+
              j] <- modX$coef[(i1 + i2 * (m + p) + i3 * (m +  p + 1)+i5*(m+p+1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                (m + p) + i3 * (m + p + 1) + i5*(m+p+1) + (j - 1) * (m + p + 1) + m + p + 1)]
      }
    }
    if(i4!=0){ # NLTD
      for (j in 1:i4){
        MAT[2:(m + p + 1), i1 + i2 + i3 + i5 + i8+j] <- modNL$coef[(i1 +
                                                                      (i2+i8) * (m + p) + (i3+i5) * (m + p + 1)+ (j - 1) * (m + p) + 1):(i1 + (i2+i8) * (m + p)
                                                                                                                                         + (i3+i5) * (m + p + 1) +  (j - 1) * (m + p) + m + p)]
        MAT[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+
              j] <- modTD$coef[(i1 + i2 * (m + p) + (i3+i5+i8) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                             (m + p) + (i3+i5+i8) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)]

      }
    }
    if(i6!=0){ # NLTEL
      for(j in 1:i6){
        MAT[2:(m + p + 1), i1 + i2 + i3 + i5 + i8+i4+j] <- modNL$coef[(i1 +
                                                                         (i2+i8+i4) * (m + p) + (i3+i5) * (m + p + 1)+ (j - 1) * (m + p) + 1):(i1 + (i2+i8+i4) * (m + p)
                                                                                                                                               + (i3+i5) * (m + p + 1) +  (j - 1) * (m + p) + m + p)]
        MAT[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+i4+
              j] <- modX$coef[(i1 + i2 * (m + p) + (i3+i5+i8) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                            (m + p) + (i3+i5+i8) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)]
      }
    }
    if(i7!=0){ #TDTEL
      for(j in 1:i7){
        MAT[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5 + i8+i6+i4+j] <- modTD$coef[(i1 +
                                                                                            (i2+i8) * (m + p) + (i3+i5+i4) * (m + p + 1)+ (j - 1) * (m + p+1) + 1):(i1 + (i2+i8) * (m + p)
                                                                                                                                                                    + (i3+i5+i4) * (m + p + 1) +  (j - 1) * (m + p+1) + m + p+1)]
        MAT[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+i6+i4+
              j] <- modX$coef[(i1 + i2 * (m + p) + (i3+i5+i8+i6) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                               (m + p) + (i3+i5+i8+i6) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)]

      }
    }



    MATse <- matrix(ncol = V, nrow = 3*(1 + m + p) ) # store the estimated SE
    if (i1 != 0) {
      for (j in 1:i1) {
        MATse[1, j] <- sqrt(diag(modX$var)[j])
      }
    }
    if (i2 != 0) {
      for (j in 1:i2) {
        MATse[2:(m + p + 1), i1 + j] <- sqrt(diag(modX$var)[(i1 + (j -
                                                                     1) * (m + p) + 1):(i1 + (j - 1) * (m + p) + m +p)])
      }
    }
    if (i3 != 0) {
      for (j in 1:i3) {
        MATse[(m + p + 2):(2 * m + 2 * p + 2), i1 + i2 + j] <- sqrt(diag(modX$var)[(i1 +
                                                                                      i2 * (m + p) + (j - 1) * (m + p + 1) + 1):(i1 +
                                                                                                                                   i2 * (m + p) + (j - 1) * (m + p + 1) + m + p + 1)])
      }
    }

    if (i5 != 0) { # only TEL no NL no TD
      for (j in 1:i5) { # row m+p+2 to 2*(m+p+1) stores TD effects
        MATse[(2*m + 2*p + 3):(3 * m + 3 * p + 3), i1 + i2 + i3+j] <-sqrt(diag(modX$var)[(i1 +
                                                                                            i2 * (m + p) + i3 * (m + p + 1) +(j - 1) * (m + p + 1)+ 1):(i1 +  i2 * (m + p) + i3 * (m + p + 1)+(j - 1) * (m + p + 1) + m + p + 1)])
      }
    }
    if(i8!=0){ #NL+TD+TEL
      for (j in 1:i8) {
        MATse[2:(m + p + 1), i1 + i2 + i3+i5 + j] <-sqrt(diag(modNL$var)[(i1 +
                                                                            i2 * (m + p) + i3 * (m + p + 1)+i5*(m+p+1) + (j - 1) * (m + p) + 1):(i1 + i2 * (m + p) + i3 * (m + p + 1) +
                                                                                                                                                   i5*(m+p+1)+(j - 1) * (m + p) + m + p)])
        MATse[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5+
                j] <-sqrt(diag(modTD$var)[(i1 + i2 * (m + p) + i3 * (m +  p + 1)+i5*(m+p+1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                            (m + p) + i3 * (m + p + 1) + i5*(m+p+1) + (j - 1) * (m + p + 1) + m + p + 1)])
        MATse[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+
                j] <- sqrt(diag(modX$var)[(i1 + i2 * (m + p) + i3 * (m +  p + 1)+i5*(m+p+1) + (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                            (m + p) + i3 * (m + p + 1) + i5*(m+p+1) + (j - 1) * (m + p + 1) + m + p + 1)])
      }
    }
    if(i4!=0){ # NLTD
      for (j in 1:i4){
        MATse[2:(m + p + 1), i1 + i2 + i3 + i5 + i8+j] <- sqrt(diag(modNL$var)[(i1 +
                                                                                  (i2+i8) * (m + p) + (i3+i5) * (m + p + 1)+ (j - 1) * (m + p) + 1):(i1 + (i2+i8) * (m + p)
                                                                                                                                                     + (i3+i5) * (m + p + 1) +  (j - 1) * (m + p) + m + p)])
        MATse[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+
                j] <- sqrt(diag(modTD$var)[(i1 + i2 * (m + p) + (i3+i5+i8) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                         (m + p) + (i3+i5+i8) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)])
      }
    }
    if(i6!=0){ # NLTEL
      for(j in 1:i6){
        MATse[2:(m + p + 1), i1 + i2 + i3 + i5 + i8+i4+j] <-sqrt(diag(modNL$var)[(i1 +
                                                                                    (i2+i8+i4) * (m + p) + (i3+i5) * (m + p + 1)+ (j - 1) * (m + p) + 1):(i1 + (i2+i8+i4) * (m + p)
                                                                                                                                                          + (i3+i5) * (m + p + 1) +  (j - 1) * (m + p) + m + p)])
        MATse[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+i4+
                j] <-sqrt(diag(modX$var)[(i1 + i2 * (m + p) + (i3+i5+i8) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                       (m + p) + (i3+i5+i8) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)])
      }
    }
    if(i7!=0){ #TDTEL
      for(j in 1:i7){
        MATse[(m + p + 2):(2 * (m + p + 1)), i1 + i2 + i3 + i5 + i8+i4+i6+j] <-sqrt(diag(modTD$var)[(i1 +
                                                                                                       (i2+i8) * (m + p) + (i3+i5+i4) * (m + p + 1)+ (j - 1) * (m + p+1) + 1):(i1 + (i2+i8) * (m + p)
                                                                                                                                                                               + (i3+i5+i4) * (m + p + 1) +  (j - 1) * (m + p+1) + m + p+1)])
        MATse[(2*m + 2*p + 3):(3 * (m + p + 1)), i1 + i2 + i3 + i5+ i8+i6+i4+
                j] <-sqrt(diag(modX$var)[(i1 + i2 * (m + p) + (i3+i5+i8+i6) * (m +  p + 1)+ (j - 1) * (m + p + 1) + 1):(i1 + i2 *
                                                                                                                          (m + p) + (i3+i5+i8+i6) * (m + p + 1) + (j - 1) * (m + p + 1) + m + p + 1)])

      }
    }


    var_order <- c(variablesNEW[NL == 0 & TD == 0 &TEL==0], variablesNEW[NL ==
                                                                           1 & TD == 0 & TEL==0], variablesNEW[NL == 0 & TD == 1&TEL==0],variablesNEW[NL == 0 & TD == 0&TEL==1],
                   variablesNEW[NL == 1 & TD == 1&TEL==1],variablesNEW[NL ==1 & TD == 1 & TEL==0],
                   variablesNEW[NL == 1 & TD == 0&TEL==1], variablesNEW[NL == 0 & TD == 1&TEL==1])


    coefficients<-MAT[1,match(variablesNEW,var_order)]
    names(coefficients)<-variables

    se_coef<-MATse[1,match(variablesNEW,var_order)]

    coefficients_splines_NL<-as.matrix(MAT[2:(m+p+1),match(variablesNEW,var_order)])
    coefficients_splines_NL<-rbind(rep(0,V),coefficients_splines_NL)
    coefficients_splines_NL[1,(NL==0)]<-NA
    colnames(coefficients_splines_NL)<-variables

    coefficients_splines_TD<-as.matrix(MAT[(m+p+2):(2*(m+p+1)),match(variablesNEW,var_order)])
    colnames(coefficients_splines_TD)<-variables

    coefficients_splines_TEL<-as.matrix(MAT[(2*(m+p)+3):(3*(m+p+1)),match(variablesNEW,var_order)])
    colnames(coefficients_splines_TEL)<-variables


    knots_covariates<-knotsNEW[1:V,]
    if (V>1) {rownames(knots_covariates)<-variables}
    if (V>1) {knots_covariates[(NL==0),]<-rep(NA,p + 1 + m + p + 1)} else { if (NL==0) {knots_covariates<-rep(NA,p + 1 + m + p + 1)} }
    knots_time<-knotsNEW[V+1,]
    if(nTEL>1){
      knots_TEL<-knotsNEW[(V+2):(V+1+nTEL),]
      rownames(knots_TEL)<-variables[TEL==1]
    }else{
      knots_TEL<-knotsNEW[(V+2),]
      names(knots_TEL)<-variables[TEL==1]
    }

    nEvents<-sum(data[,TypeNEW[3]]==1)
    print(vrais)
    rm(data)
    gc()
    list(Partial_Log_Likelihood = vrais[length(vrais)],
         log_likelihood_history = vrais,Number_of_parameters = nbpara +
           nbonlyNL * (m + p) + (nbonlyTD+nbonlyTEL+2*nbTDTEL) * (1 + m + p) + (nbNLTD+nbNLTEL) *
           (m + p + m + p + 1)+nbNLTDTEL*(3*(m+p+1)-1), Number_events=nEvents, Number_knots = m, Degree_of_splines = p,
         knots_covariates = knots_covariates,
         knots_time = knots_time,
         knots_TEL = knots_TEL,
         coefficients = coefficients, Standard_Error=se_coef,
         coefficients_splines_NL = coefficients_splines_NL,
         coefficients_splines_TEL = coefficients_splines_TEL,
         coefficients_splines_TD = coefficients_splines_TD,variables=variables)
  }

}
