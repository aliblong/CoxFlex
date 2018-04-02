DvlpMatrix<-function(data, listeT, ncol, TypeNEW) {
  data <- matrix(data, ncol = ncol)
  if (max(data[, TypeNEW[2]]) < min(listeT[listeT != 0])) {
    XX <- data
  }
  else {
    aindex <- rep(0, (sum(listeT <= max(data[, TypeNEW[2]])) -
                        1))

    for (i in 1:(sum(listeT <= max(data[, TypeNEW[2]])) -
                 1)) {
      for (j in 1:(dim(data)[1])) {

        if (as.numeric(data[j, TypeNEW[1]]) < as.numeric(listeT[1 + i]) & as.numeric(data[j, TypeNEW[2]]) >= as.numeric(listeT[1 + i]))
          aindex[i] <- j
      }
    }
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
