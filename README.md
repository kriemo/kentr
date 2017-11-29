
kentr 
====
![Build Status](https://travis-ci.org/kriemo/kentr.svg?branch=master)

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

``` r
library(kentr)
fa_path <- system.file("extdata", "test.fasta", package = "kentr")
df <- data.frame(chrom = "chr1",
                 start = c(20000, 25000),
                 end = c(20025, 25025))

df
#>   chrom start   end
#> 1  chr1 20000 20025
#> 2  chr1 25000 25025

seqs <- get_sequences(df, fa_path)

seq1 <- seqs[1, "seq"]
seq2 <- seqs[2, "seq"]
get_hamming(seq1, seq2)
#> [1] 14
```

To count kmers in sequences:
``` r
library(kentr)
library(tidyverse)

fa_path <- system.file("extdata", "test.fasta", package = "kentr")
df <- data.frame(chrom = "chr1",
                 start = c(20000, 25000),
                 end = c(20500, 25500))

df
#>   chrom start   end
#> 1  chr1 20000 20500
#> 2  chr1 25000 25500

seqs <- get_sequences(df, fa_path)

seqs
#>   chrom start   end           header
#> 1  chr1 20000 20500 chr1:20000-20500
#> 2  chr1 25000 25500 chr1:25000-25500
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    seq
#> 1 cctggtgctcccacaaaggagaagggctgatcactcaaagttgcgaacaccaagctcaacaatgagccctggaaaatttctggaatggattattaaacagagagtctgtaagcacttagaaaaggccgcggtgagtcccaggggccagcactgctcgaaatgtacagcatttctctttgtaacaggattattagcctgctgtgcccggggaaaacatgcagcacagtgcatctcgagtcagcaggattttgacggcttctaacaaaatcttgtagacaagatggagctatgggggttggaggagagaacatataggaaaaatcagagccaaatgaaccacagccccaaagggcacagttgaacaatggactgattccagccttgcacggagggatctggcagagtCCATCCAGTTCATTCAACACCTGGTTAGAAAACTGGGGCCAGCACACAGGGGAAGGGTAAGCTGGTTTCATGATCGAATCAAGGCTCAGACAATT
#> 2 GCTTCAGCCTGCACAGATAGGGGAGTAGGGGACAGAGCATTTGCTGAGAGGCCAGGAGCGCATAGATGGGACTCTGCTGATGCCTGCTGAGTGAATGAGGGAAAGGGCAGGGCCCGGGACTGGGGAATCTGTAGGGTCAATGGAGGAGTTCAGAGAAGGTGCAACATTTCTGACCCCCTACAAGGTGCTTGCTACCTGCCAGGCACCCTTTCCATACCTTGTCTCAGTTCAGCTCCCCACCTTGGATAAACAAGAAACCTTGGTTGCAGAGGAAAAAAGAGGCTGGAAACAAAGGGGTAGAAATGGGGTAGCAGGGGAGATTGCCTGATCAACTGCCAAATGGTACACAGTTCTGGAAAAGCACAAAAAATGTGCACACACGGGTTCTTCCCACTTTAACCCCTGAGGAATCTGAGGCCTGCTCCTGAAACAGACTGGGCAGTGGCTAGTGACTCTAGGTATAGGAGTATCCAGCCCTGCTCACCCAGGCTAGAGCTTAG

# returns list with a dataframe of kmers and counts for each sequence
get_kmers(seqs$seq, n = 2)
#> [[1]]
#>    kmer counts
#> 1    AA     51
#> 2    AC     28
#> 3    AG     47
#> 4    AT     29
#> 5    CA     50
#> 6    CC     24
#> 7    CG      9
#> 8    CT     27
#> 9    GA     39
#> 10   GC     31
#> 11   GG     44
#> 12   GT     19
#> 13   TA     15
#> 14   TC     26
#> 15   TG     33
#> 16   TT     27
#> 
#> [[2]]
#>    kmer counts
#> 1    AA     39
#> 2    AC     27
#> 3    AG     51
#> 4    AT     19
#> 5    CA     38
#> 6    CC     34
#> 7    CG      3
#> 8    CT     40
#> 9    GA     40
#> 10   GC     34
#> 11   GG     54
#> 12   GT     20
#> 13   TA     19
#> 14   TC     20
#> 15   TG     40
#> 16   TT     21

# use with dplyr tibbles 

seqs <- as_data_frame(seqs)
seqs
#> # A tibble: 2 x 5
#>    chrom start   end           header
#>   <fctr> <dbl> <dbl>            <chr>
#> 1   chr1 20000 20500 chr1:20000-20500
#> 2   chr1 25000 25500 chr1:25000-25500
#> # ... with 1 more variables: seq <chr>

kmers <- mutate(seqs, 
                kmers = get_kmers(seq)) %>% 
                select(header, kmers) 
                
kmers
#> # A tibble: 2 x 2
#>             header                 kmers
#>              <chr>                <list>
#> 1 chr1:20000-20500 <data.frame [16 x 2]>
#> 2 chr1:25000-25500 <data.frame [16 x 2]>

unnest(kmers)
#> # A tibble: 32 x 3
#>              header   kmer counts
#>               <chr> <fctr>  <int>
#>  1 chr1:20000-20500     AA     51
#>  2 chr1:20000-20500     AC     28
#>  3 chr1:20000-20500     AG     47
#>  4 chr1:20000-20500     AT     29
#>  5 chr1:20000-20500     CA     50
#>  6 chr1:20000-20500     CC     24
#>  7 chr1:20000-20500     CG      9
#>  8 chr1:20000-20500     CT     27
#>  9 chr1:20000-20500     GA     39
#> 10 chr1:20000-20500     GC     31
#> # ... with 22 more rows
```
