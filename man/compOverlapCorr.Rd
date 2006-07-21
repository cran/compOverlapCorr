\name{compOverlapCorr}
\alias{compOverlapCorr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Comparing overlapping correlation coefficients }
\description{
  This package contains function to test the difference between two overlapping (in the sense of having a variable in common) correlation coefficients, using a Z-test as described by Meng, Rosenthal, and Rubin (Meng et al, 1992).
}
\usage{
compOverlapCorr(N, r13, r23, r12)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{N}{ the sample size }
  \item{r13}{ the correlation coefficient between the two variables, X1 and Y }
  \item{r23}{ the correlation coefficient between the two variables, X2 and Y }
  \item{r12}{ the correlation coefficient between the two variables, X1 and X2 }
}
\value{
  \item{comp1 }{the Z-score generated using the Meng-Rosenthal-Rubin method}
  \item{comp2 }{the P-value, two-tailed, for the difference between r13 and r23}
  
}
\references{ Meng X-L, Rosenthal R, Rubin DB. Comparing correlated correlation coefficients. Psychological bulletin 1992;111(1):172-175.}
\author{ Ka-Loh Li <ka-loh.li@radiology.ucsf.edu>, and Xiaoping Zhu <Xiaoping.Zhu@ucsf.edu> }
\examples{resul <- compOverlapCorr(15, 0.63, -0.03, -0.19)


## The function is currently defined as
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
}
\keyword{ methods }% at least one, from doc/KEYWORDS

