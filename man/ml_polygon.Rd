% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ml_polygon.R
\name{ml_polygon}
\alias{ml_polygon}
\title{ml_polygon}
\usage{
ml_polygon(
  bufferPoints.sf,
  polygons.sf,
  polygons.id = "SPAID",
  weightsCol = "layer"
)
}
\arguments{
\item{bufferPoints.sf}{Buffer points of a given displaced location containing weights}

\item{polygons.sf}{Polygons to join the displaced points to. For now, the function assumes there is a column with unique ID for each polygon called "SPAID"}

\item{polygons.id}{Name of column of polygons.id sf object that contains a unique ID, default is "SPAID".}

\item{weightsCol}{Name of the column in the bufferPoints.sf object that contains the weights to use, defaults to "layer".}
}
\value{
returns a dataframe with one row containing the data from most likely health facility joined, and its estimated probability.
}
\description{
Determines the maximum likely polygon for a given point
}
\details{
Replaces the weightedClose function
}
\examples{

# coming soon!


}
\author{
J.W. Rozelle
}
