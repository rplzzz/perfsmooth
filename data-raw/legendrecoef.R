#### Coefficients for the Legendre polynomials through degree 10
### NB: This could easily be extended to higher degree, but one needs to be
###     careful about loss of precision in the coefficient formula for higher degree

dmax <- 10
legendrecoef <- matrix(0, nrow=dmax+1, ncol=dmax+1)

for(n in seq(0,dmax)) {
  a <- rep(0, dmax+1)
  for(k in seq(0,n)) {
    ik <- k+1
    a[ik] <- 2^n * choose(n,k) * choose((n+k-1)/2, n)
  }
  legendrecoef[,n+1] <- a
}

colnames(legendrecoef) <- paste0('P', seq(0,dmax))
row.names(legendrecoef) <- paste0('x^{',seq(0,dmax),'}')

usethis::use_data(legendrecoef, overwrite = TRUE)
