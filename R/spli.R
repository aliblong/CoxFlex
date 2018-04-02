spli<-function(x, j, p, knots) {
  if (p == 0) {
    b <- ifelse(x >= knots[j] & x < knots[j + 1], 1,
                0)
    return(b)
  }
  else {
    a1 <- ifelse(rep(knots[j] != knots[j + p], length(x)),
                 (x - knots[j])/(knots[j + p] - knots[j]), 0)
    a2 <- ifelse(rep(knots[j + p + 1] != knots[j + 1],
                     length(x)), (knots[j + p + 1] - x)/(knots[j +
                                                                 p + 1] - knots[j + 1]), 0)
    return(a1 * spli(x, j, p - 1, knots) + a2 * spli(x,
                                                     j + 1, p - 1, knots))
  }
}
