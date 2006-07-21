"compOverlapCorr" <-
function(N, r13, r23, r12){
# Meng-Rosenthal-Rubin method 
# comparing two correlated correlation coefficients (r13, r23)
# N is the number of subjects
# r13 is the correlation coefficient between variables X1 and Y
# r23 is the correlation coefficient between variables X2 and Y 
# r12 is the correlation coefficient between variables X1 and X2 
# return difference and two-tailed p-value as list(diff, pval) 

# Fisher Z-transform 
  zf1 <- 0.5*log((1 + r13)/(1 - r13)) 
  zf2 <- 0.5*log((1 + r23)/(1 - r23))

# difference
  r2mean <- (r13^2 + r23^2)/2.0
  f <- (1-r12)/(2*(1-r2mean))
  if(f > 1.0) 
        {
        f <- 1.0
        }
  h <- (1-f*r2mean)/(1-r2mean)
  dz <- (zf1 - zf2)*sqrt((N-3)/(2*(1-r12)*h))

# two-tailed p-value 
  pv <- 2*(1 - pnorm(abs(dz)))
  #pv <- 1 - pnorm(abs(dz))
  return(as.numeric(list(diff=dz, pval=pv))) 
}

