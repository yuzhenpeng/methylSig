% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methylSig_data_doc.R
\docType{data}
\name{sample_data}
\alias{sample_data}
\title{Sample data}
\format{A mixture of BSseq and GRanges class objects. See documentation for test_data.}
\description{
The \code{sample_data} object contains the following items:
\describe{
  \item{meth}{ A \code{BSseq-class} object containing 6 samples total, with three in each group. Genome is hg19. }
  \item{tiled_meth}{ A tiled version of the \code{BSseq-class} object called \code{meth}. Tiles are 1000bp. Genome is hg19. }
  \item{msig_cpgs}{ A \code{GRanges-class} object containing the results of \code{methylSigCalc} on \code{data}. }
  \item{msig_tiles}{ A \code{GRanges-class} object containing the results of \code{methylSigCalc} on \code{tiled_meth}. }
  \item{tfbs}{ A \code{GRanges-class} object representing transcription factor binding sites. For use in \code{methylSig.tfbsEnrichTest}. Genome is hg19. }
}
}
\keyword{datasets}
