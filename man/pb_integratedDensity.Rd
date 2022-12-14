% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pb_integratedDensity.R
\name{pb_integratedDensity}
\alias{pb_integratedDensity}
\title{pb_integratedDensity}
\usage{
pb_integratedDensity(
  cellMeters = 50,
  boundaries = c(5000, 10000),
  weights = c(0.99, 0.01)
)
}
\arguments{
\item{cellMeters}{Resolution of the probability buffer. Higher resolutions will result in more accuracy, but will take more computational resources to use. Default is 50 meters.}

\item{boundaries}{A numeric vector of the radius limits. Defaults to 5,000 and 10,000, which is the displacement for rural clusters in the DHS (99% displaced between 0 and 5 km, 1% displaced between 0 and 10 km)}

\item{weights}{A numeric vector for weights of the boundaries set in the boundaries argument, if using more than one boundary. Defaults to c(0.99, 0.01) based on the rural displacement probability}
}
\value{
Returns a list with a matrix of the weights from the pdf, the input resolution, and the input boundaries
}
\description{
This function integrates the probability density function for raster based on provided resolution, radii, and weights. Generally this function is not useful by itself, but is used as a helper function.
}
\examples{

# create the probability density matrix
ruralPB <- pb_integratedDensity()

ruralPB.matrix <- ruralPB[["weightedCircle"]]


}
\author{
J.W. Rozelle
}
