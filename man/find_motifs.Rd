% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_motifs.R
\name{find_motifs}
\alias{find_motifs}
\title{Find positions of query seq matches in seq column}
\usage{
find_motifs(df, query_seq)
}
\arguments{
\item{df}{data frame containing a header and sequence column}

\item{query_seq}{sequence to identify in df$seq column, can contain regex}
}
\value{
A data.frame with chrom, start, end, and header columns.
Start and end indicate start
point and end point of query_seq match.
}
\description{
Find positions of query seq matches in seq column
}
\examples{
fa_path <- system.file("extdata", "test.fasta", package = "kentr")
query_seq <- "AATAAA[GTC]"

seqs <- read_genome_seq(fa_path)
find_motifs(seqs, query_seq)

}
