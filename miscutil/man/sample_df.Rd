% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sample_df.R
\name{sample_df}
\alias{sample_df}
\title{Random sampling of a data frame}
\usage{
sample_df(x, size.row = ifelse(unit == "total", nrow(x), 100),
  size.col = ifelse(unit == "total", ncol(x), 100), unit = "total",
  replace.row = F, replace.col = F)
}
\arguments{
\item{x}{data.frame or similar object.}

\item{size.row}{numeric. number of rows if \code{unit == 'total'}, otherwise percentage of rows.}

\item{size.col}{numeric. number of columns if \code{unit == 'total'}, otherwise percentage of columns.}

\item{unit}{either 'total' or 'percent'.}

\item{replace.row}{logical. should rows be sampled with replacement?}

\item{replace.col}{logical. should columns be sampled with replacement?}
}
\value{
a data frame with sampled rows and/or columns.
}
\description{
\code{sample_df} takes a specified number of rows and/or columns from \code{x} either with or without replacement.
}
\examples{
x = matrix(1:500, nrow = 50, dimnames = list(paste0('r',1:50),paste0('c',1:10)))
x1 = sample_df(x, size.row = 4)
x2 = sample_df(x, size.row = 4, unit = 'percent')
x3 = sample_df(x, size.row = 4, size.col = 5, unit = 'total')
}
