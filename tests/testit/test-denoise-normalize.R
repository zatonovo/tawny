assert('normalize returns', {
  x <- matrix(rnorm(80, sd=4), ncol=4)
  y <- normalize(x)
  v <- apply(y, 2, var)

  all.equal(v, rep(1,4), tolerance=0.0001)
})
