
# kentr

![R-CMD-check](https://github.com/kriemo/kentr/workflows/R-CMD-check/badge.svg)

A R package with a collection of functions for working with sequence
data in R. Includes the
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
library(tidyverse, warn.conflicts = FALSE)

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

seqs <- as_tibble(seqs)
seqs
#> # A tibble: 2 x 5
#>   chrom start   end header       seq                                            
#>   <fct> <dbl> <dbl> <chr>        <chr>                                          
#> 1 chr1  20000 20500 chr1:20000-… cctggtgctcccacaaaggagaagggctgatcactcaaagttgcga…
#> 2 chr1  25000 25500 chr1:25000-… GCTTCAGCCTGCACAGATAGGGGAGTAGGGGACAGAGCATTTGCTG…

kmers <- mutate(seqs, 
                kmers = get_kmers(seq)) %>% 
                select(header, kmers) 
                
kmers
#> # A tibble: 2 x 2
#>   header           kmers            
#>   <chr>            <list>           
#> 1 chr1:20000-20500 <df[,2] [16 × 2]>
#> 2 chr1:25000-25500 <df[,2] [16 × 2]>

unnest(kmers)
#> Warning: `cols` is now required.
#> Please use `cols = c(kmers)`
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
#> # … with 22 more rows
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

### Extract bam alignments as a data.frame

``` r

bam_file <- system.file("extdata",
                        "small_sorted.bam", 
                        package = "kentr")
res <- bam_to_df(bam_file)
head(res)
#>   chrom   start     end                                            name strand
#> 1  chr1 3015399 3015499 HWI-ST1133R:8:2304:12130:31201#GAACCCATTCTTTCCC      -
#> 2  chr1 3015400 3015500 HWI-ST1133R:8:1313:13728:38795#GAACCCATTCTTTCCC      +
#> 3  chr1 3019635 3019735 HWI-ST1133R:8:1209:10428:12865#GAACCCATTCTTTCCC      +
#> 4  chr1 3019635 3019735  HWI-ST1133R:8:2315:2419:54036#GAACCCATTCTTTCCC      -
#> 5  chr1 3043511 3043611  HWI-ST1133R:8:1107:9218:64408#GAACCCATTCTTTCCC      -
#> 6  chr1 3053963 3054063  HWI-ST1133R:8:2212:2473:90987#GAACCCATTCTTTCCC      +

res <- bam_to_df(bam_file,
                 region = "chr1:3e6-3.1e6",
                 tags = c("XS:A", "AS:i"))
res
#>    chrom   start     end                                            name strand
#> 1   chr1 3015399 3015499 HWI-ST1133R:8:2304:12130:31201#GAACCCATTCTTTCCC      -
#> 2   chr1 3015400 3015500 HWI-ST1133R:8:1313:13728:38795#GAACCCATTCTTTCCC      +
#> 3   chr1 3019635 3019735 HWI-ST1133R:8:1209:10428:12865#GAACCCATTCTTTCCC      +
#> 4   chr1 3019635 3019735  HWI-ST1133R:8:2315:2419:54036#GAACCCATTCTTTCCC      -
#> 5   chr1 3043511 3043611  HWI-ST1133R:8:1107:9218:64408#GAACCCATTCTTTCCC      -
#> 6   chr1 3053963 3054063  HWI-ST1133R:8:2212:2473:90987#GAACCCATTCTTTCCC      +
#> 7   chr1 3054014 3054114  HWI-ST1133R:8:2212:2473:90987#GAACCCATTCTTTCCC      -
#> 8   chr1 3058369 3058469   HWI-ST1133R:8:2309:7368:5858#GAACCCATTCTTTCCC      -
#> 9   chr1 3058645 3058745 HWI-ST1133R:8:1111:14310:29996#GAACCCATTCTTTCCC      -
#> 10  chr1 3058647 3058746 HWI-ST1133R:8:2302:13219:91951#GAACCCATTCTTTCCC      +
#> 11  chr1 3059365 3059465  HWI-ST1133R:8:1212:8780:42946#GAACCCATTCTTTCCC      +
#> 12  chr1 3068842 3068939 HWI-ST1133R:8:2215:20128:28486#GAACCCATTCTTTCCC      +
#> 13  chr1 3068851 3068939 HWI-ST1133R:8:2215:20128:28486#GAACCCATTCTTTCCC      -
#> 14  chr1 3077199 3077299  HWI-ST1133R:8:1303:9194:69620#GAACCCATTCTTTCCC      -
#> 15  chr1 3090805 3090905  HWI-ST1133R:8:1216:5060:21458#GAACCCATTCTTTCCC      +
#> 16  chr1 3090822 3090922  HWI-ST1133R:8:1216:5060:21458#GAACCCATTCTTTCCC      -
#>    XS  AS
#> 1  45   0
#> 2  45   0
#> 3  45 -18
#> 4  45 -18
#> 5  43   0
#> 6  45   0
#> 7  45   0
#> 8  43 -12
#> 9  45  -6
#> 10 45 -13
#> 11 45  -6
#> 12 43 -15
#> 13 43 -18
#> 14 43   0
#> 15 45   0
#> 16 45   0
```

### Find motif matches genome-wide

``` r
fa_path <- system.file("extdata", "test.fasta", package = "kentr")
query_seq <- "AATAAA[GTC]"

matches <- read_genome_seq(fa_path) %>% 
  find_motifs(query_seq) %>% 
  head()

matches
#>   chrom start   end        header
#> 1  chr1 11544 11551 chr1:0-100000
#> 2  chr1 33977 33984 chr1:0-100000
#> 3  chr1 34608 34615 chr1:0-100000
#> 4  chr1 35489 35496 chr1:0-100000
#> 5  chr1 35497 35504 chr1:0-100000
#> 6  chr1 37324 37331 chr1:0-100000

get_sequences(matches, fa_path)
#>   chrom start   end        header          header1     seq
#> 1  chr1 11544 11551 chr1:0-100000 chr1:11544-11551 aataaat
#> 2  chr1 33977 33984 chr1:0-100000 chr1:33977-33984 aataaac
#> 3  chr1 34608 34615 chr1:0-100000 chr1:34608-34615 aataaag
#> 4  chr1 35489 35496 chr1:0-100000 chr1:35489-35496 aataaat
#> 5  chr1 35497 35504 chr1:0-100000 chr1:35497-35504 aaTAAAT
#> 6  chr1 37324 37331 chr1:0-100000 chr1:37324-37331 aataaat
```
