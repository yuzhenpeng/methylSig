% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read.R
\name{methylSigReadData}
\alias{methylSigReadData}
\title{Read methylation score files to make a 'BSseq' object.}
\usage{
methylSigReadData(fileList, pData, assembly = NA, destranded = TRUE,
  maxCount = 500, minCount = 10, filterSNPs = FALSE, num.cores = 1,
  fileType = c("cov", "cytosineReport"), verbose = TRUE)
}
\arguments{
\item{fileList}{Files to be read. These can be \code{cov} or \code{cytosine_reports} from the Bismark Methylation Extractor. See \code{fileType} for details.}

\item{pData}{A \code{data.frame} containing phenotype information for the samples in \code{fileList}. The \code{row.names} attribute of the \code{data.frame} should match the \code{Sample_Names}. See example below.}

\item{assembly}{The genome assembly used for alignment. e.g. \code{hg19}, \code{mm10}, etc.}

\item{destranded}{A logical value indicating whether to destrand the reverse to forward strand. If TRUE, the reads from both will be combined. Default is TRUE.}

\item{maxCount}{A number indicating the maximum coverage count to be included.}

\item{minCount}{A number indicating the minimum coverage count to be included.}

\item{filterSNPs}{A logical value indicating whether or not to filter out C > T SNPs based on the 1000 Genomes Project. NOTE: Only supported when \code{assembly = 'hg19'}.}

\item{num.cores}{Number of cores to be used in reading files. Default is 1.}

\item{fileType}{The format of the input file. Either \code{cov} or \code{cytosineReport}. One of the outputs of the Bismark Methylation Extractor.}

\item{verbose}{A logical value indicating whether \code{bsseq::read.bismark} shoud print progress. Default TRUE.}
}
\value{
A \code{BSseq-class} object.
}
\description{
This function reads files created by the Bismark Methylation Extractor, and outputs a \code{BSseq} object.
}
\examples{
files = c(
    system.file('extdata', 'MDAMB_231_1DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_1DS.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_2DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_2DS.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_3DR.txt.gz', package='methylSig'),
    system.file('extdata', 'MDAMB_231_3DS.txt.gz', package='methylSig'))

sample.ids = basename(files)
sample.ids = gsub('.txt.gz', '', sample.ids)

pData = data.frame(
    Sample_Names = sample.ids,
    DR_vs_DS = relevel(factor(c('DR','DS','DR','DS','DR','DS')), ref = 'DS'),
    row.names = sample.ids,
    stringsAsFactors = FALSE)

meth = methylSigReadData(
    fileList = files,
    pData = pData,
    assembly = 'hg19',
    destranded = TRUE,
    maxCount = 500,
    minCount = 10,
    filterSNPs = TRUE,
    num.cores = 1,
    fileType = 'cytosineReport')

}
\seealso{
\code{\link{methylSigCalc}}
}
