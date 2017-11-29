
kentr
====

An R package with a collection of functions for working with sequence data in R. 
Wraps the [`htslib`](http://www.htslib.org/) C library for querying indexed fasta files.

```r
#install.packages('devtools')
devtools::install_github('kriemo/kentr')
```

Basic Usage
===========

To extract sequences from indexed fasta:
``` r
library(kentr)
# generate data.frame (or tibble) with bed coordinates

df <- data.frame(chrom = "chr1",
                 start = 20000,
                 end = 20025)

df
#>   chrom start   end
#> 1  chr1 20000 20025

# path to fasta file (with fai samtools faidx index)
fa_path <- system.file("extdata", "test.fasta", package = "kentr")

get_sequences(df, fa_path)
#>   chrom start   end           header                       seq
#> 1  chr1 20000 20025 chr1:20000-20025 cctggtgctcccacaaaggagaagg
```

To compute hamming distance between sequences:

```r
fa_path <- system.file("extdata", "test.fasta", package = "kentr")
df <- data.frame(chrom = "chr1",
                 start = c(20000, 25000),
                 end = c(20025, 25025))

df

seqs <- get_sequences(df, fa_path)

seq1 <- seqs[1, "seq"]
seq2 <- seqs[2, "seq"]
get_hamming(seq1, seq2)

```

To count kmers in sequences
