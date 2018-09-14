
# kentr

[![Build
Status](https://travis-ci.org/kriemo/kentr.svg?branch=master)](https://travis-ci.org/kriemo/kentr)

A R package with a collection of functions for working with sequence
data in R. Wraps the [`htslib`](http://www.htslib.org/) C library and
the
[`ssw`](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
Smith-Waterman alignment C/C++ library.

``` r
#install.packages('devtools')
devtools::install_github('kriemo/kentr')
```

## Basic Usage

### Extract sequences from fasta:

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

### Compute Hamming distances:

``` r
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

### Count kmers in sequences

``` r
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
#>   chrom start   end header       seq                                      
#>   <fct> <dbl> <dbl> <chr>        <chr>                                    
#> 1 chr1  20000 20500 chr1:20000-… cctggtgctcccacaaaggagaagggctgatcactcaaag…
#> 2 chr1  25000 25500 chr1:25000-… GCTTCAGCCTGCACAGATAGGGGAGTAGGGGACAGAGCAT…

kmers <- mutate(seqs, 
                kmers = get_kmers(seq)) %>% 
                select(header, kmers) 
                
kmers
#> # A tibble: 2 x 2
#>   header           kmers                
#>   <chr>            <list>               
#> 1 chr1:20000-20500 <data.frame [16 × 2]>
#> 2 chr1:25000-25500 <data.frame [16 × 2]>

unnest(kmers)
#> # A tibble: 32 x 3
#>    header           kmer  counts
#>    <chr>            <chr>  <int>
#>  1 chr1:20000-20500 AA        51
#>  2 chr1:20000-20500 AC        28
#>  3 chr1:20000-20500 AG        47
#>  4 chr1:20000-20500 AT        29
#>  5 chr1:20000-20500 CA        50
#>  6 chr1:20000-20500 CC        24
#>  7 chr1:20000-20500 CG         9
#>  8 chr1:20000-20500 CT        27
#>  9 chr1:20000-20500 GA        39
#> 10 chr1:20000-20500 GC        31
#> # ... with 22 more rows
```

### Perform Smith-Waterman alignment

Uses the
[Complete-Striped-Smith-Waterman-Library](https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library)
to perform alignment between a query sequence and vector of reference
sequences.

``` r
query_seq <- seqs$seq[1]
ref_seqs <- seqs$seq

get_sw(query_seq, ref_seqs)
#>   sw_score secondary_sw_score rstart rend qstart qend secondary_rend
#> 1     1000                498      0  499      0  499            248
#> 2      112                 49     11  482     11  492            190
#>   nmismatches
#> 1           0
#> 2         272
#>                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   cigar
#> 1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  500=
#> 2 11S4=1X1=4D4=1X4=2I2=5I2=1X2=3D4=1D2=1D1=2X3=2D3=1X2=3X3=1X1=1D1=1X3=2D1=3I2=2X3=2D2=1X2=1X2=1X1=3X2=1I2=3D2=1X1=1X1=1X1=3I1=1X3=2X1=1X2=5I1X1=1X2=1D1=1X2=1I1=1X1=1X2=3I2=1X1=2I1=1X2=1D2=1D2=1X2=1X2=1X7=6D2=5I3=1X1=1X1=5I2=1X4=2D1X4=1X2=6I2=5I1=1X1=2X1=1X3=2D2=1X1=2X4=5I2=4I1=1X3=1X1=2X2=5I4=2D2=1X6=2I1X1=1X2=1I1=3I5=2D4=5D4=6D4=1I1=1I3=3D1=1X4=7D2=3D2=4D3=3I3=1X1=1X2=1X1=2I5=1X2=2D7=3D3=1I2=1I1=1X2=2X1=4I1=4X5=2X1=1D3=1X3=7I2=3I3=3I1=1I2=1X3=1X4=3I2=1X2=1X3=1X4=1X2=5D1=1X4=3D3=1X2=1D2=2X2=3D2=3D3=1X1=1I1=1X2=2I1=1X3=1X2=4D5=7S
#>   strand
#> 1      +
#> 2      +
```
