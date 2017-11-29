#'Extract sequences from fasta using bed intervals
#'
#'@param df dataframe containing bed intervals (columm names chrom, start, end, and optionally strand)
#'@param fasta_path path to fasta file with .fai index
#'@param strand return sequences based on strand (default is TRUE if strand exists)
#'@return input dataframe is returned with additional columns `header` and `seq` containing
#'  a header for with the queried interval in samtools region format (1based) and the DNA sequence
#'@export
get_sequences <- function(df, fasta_path, strand = TRUE){

  stopifnot(all(c("chrom", "start", "end") %in% colnames(df)))

  seqs <- getSeq(df, path.expand(fasta_path))
  res <- dplyr::bind_cols(df, seqs)

  if(strand & ("strand" %in% colnames(df))){
    res <- mutate(res,
                  seq = ifelse(strand == "+",
                               seq,
                               revComp(seq)),
                  header = paste0(header, "(", strand, ")"))

  }
  res
}

#'Write fasta format data to disk
#'@param df dataframe containing header and sequence for fasta
#'@param out_path filename for output fasta records
#'@param header_col column name containing header
#'@param seq_col column name containing sequences
#'@param gz use [R.utils::gzip] to compress output and add `.gz` suffix
#'@export
write_fasta <- function(df, out_path,
                        header_col = "header",
                        seq_col = "seq",
                        gz = FALSE){

  stopifnot(all(c(header_col, seq_col) %in% colnames(df)))

  res <- paste0(">", df[[header_col]], "\n", df[[seq_col]])
  write_lines(res, path.expand(out_path))

  if(gz) R.utils::gzip(path.expand(out_path), remove = T, overwrite = T)

}
