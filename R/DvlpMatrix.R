DvlpMatrix<-function(data, listeT, ncol, TypeNEW) {
  data <- matrix(data, ncol = ncol) # data is a list of full data on each individual, first transform this list to the matrix
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])) {
    XX <- data  ## if max(stop time) is less then min(event time), then data remain the same
  }
  else { ## if max(stop time of certain individual) > min(event time)
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) -
                        1))
    # apb <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) -
    #                  1))
    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) -  ## for the length of unique event time
                 1)) {
      for (j in 1:(dim(data)[1])) { #for each row of certain individual, if start time < ith order of event time<=stop time of row j, stores
        # the largest row satisfy this condition in aindex[i]
        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1 + i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1 + i]))
          aindex[i] <- j
      }
    }   ## xx first coloum repeat ID, coloum for start day assigns c(0,unique event time-last one), colum for stop day assigns unique event time
    XX <- matrix(nrow = sum(listeT <= max(data[, TypeNEW[2]])) -
                   1, ncol = ncol)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), 1] <- rep(data[1,
                                                                      1], sum(listeT <= max(data[, TypeNEW[2]])) - 1)
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[1]] <- listeT[1:(sum(listeT <=
                                                                                      max(data[, TypeNEW[2]])) - 1)]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[2]] <- listeT[2:(sum(listeT <=
                                                                                      max(data[, TypeNEW[2]])))]
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), TypeNEW[3]] <- c(rep(0,(sum(listeT <= max(data[, TypeNEW[2]])) - 2)), data[dim(data)[1],TypeNEW[3]])
    XX[1:(sum(listeT <= max(data[, TypeNEW[2]])) - 1), -c(1,
                                                          TypeNEW[1], TypeNEW[2], TypeNEW[3])] <- as.matrix(data[aindex,
                                                                                                                 -c(1, TypeNEW[1], TypeNEW[2], TypeNEW[3])])

  }
  X <- XX
  list(X)
}
