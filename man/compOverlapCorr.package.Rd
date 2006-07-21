\name{compOverlapCorr-package}
\alias{compOverlapCorr-package}
\alias{compOverlappingCorr}
\docType{package}
\title{Comparing overlapping correlation coefficients
}
\description{This package contains function to test the difference between two overlapping (in the sense of having a variable in common) correlation coefficients, using a Z-test as described by Meng, Rosenthal, and Rubin (Meng et al, 1992).
}
\details{
\tabular{ll}{
Package: \tab compOverlapCorr\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2006-07-21\cr
License: \tab GPL version 2\cr
}
}
\author{
Ka-Loh Li <ka-loh.li@radiology.ucsf.edu>, and Xiaoping Zhu <Xiaoping.Zhu@ucsf.edu>

Maintainer: Ka-Loh Li <ka-loh.li@radiology.ucsf.edu>}
\references{
Meng X-L, Rosenthal R, Rubin DB. Comparing correlated correlation coefficients. Psychological bulletin 1992;111(1):172-175.
}
\keyword{ package }
\examples{resul <- compOverlapCorr(15, 0.63, -0.03, -0.19)}
