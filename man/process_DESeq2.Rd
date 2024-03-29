% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_DESeq2.R
\name{process_DESeq2}
\alias{process_DESeq2}
\title{Process the DAF analysis through the DESeq2 package}
\usage{
process_DESeq2(data, ...)
}
\arguments{
\item{data}{the ouput of the \code{\link{build_DAF_data}} function}

\item{...}{additionnal parameters of the method}
}
\value{
the output of the DESeq() function processed through the results() function
}
\description{
Process the DAF analysis through the DESeq2 package
}
