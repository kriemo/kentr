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

#'Read fasta format into data.frame
#'@param fasta_path filename for input fasta records
#'@importFrom readr read_file
#'@importFrom stringr str_split str_trim
#'@importFrom tidyr separate
#'@export
read_fasta <- function(fasta_path){

  fasta_lines <- readr::read_file(fasta_path)
  fasta_lines <- stringr::str_split(fasta_lines, ">")[[1]]
  fasta_lines <- fasta_lines[fasta_lines != ""]
  fasta_lines <- stringr::str_trim(fasta_lines)
  df <- data.frame(records = fasta_lines)
  df <- tidyr::separate(df, .data$records, into = c("name", "seq"), sep = "\n")
  df
}

#'Extract sequences from fasta using gtf input
#'
#'@param df dataframe containing gtf style entries, either in as bed (chrom, start, end) or (seqnames)
#'@param fasta_path path to fasta file with .fai index
#'@return data_frame containing transcript id column and sequence column
#'@export
gtf_to_seq <- function(df, fasta_path){

  df_cols <- c("chrom", "start", "end",
               "transcript_id", "exon_number", "strand")

  gtf_cols <- c("seqnames", "start", "end",
                "transcript_id", "exon_number", "strand")

  if (!all(df_cols %in% colnames(df))) {
    if (all(gtf_cols %in% colnames(df))){
      df <- dplyr::rename(df, chrom = seqnames)
      df <- dplyr::mutate(df, start = start - 1)
    } else {
      stop("unknown columns in supplied dataframe")
    }
  }

  df <- dplyr::filter(df, type == "exon")
  df <- dplyr::select(df,
                      chrom,
                      start,
                      end,
                      transcript_id,
                      exon_number,
                      strand)

  seq_df <- get_sequences(df, fasta_path)
  seq_df <- dplyr::group_by(seq_df, transcript_id)
  seq_df <- dplyr::arrange(seq_df, exon_number, .by_group = TRUE)
  res <- dplyr::summarize(seq_df, seq = stringr::str_c(seq, collapse = ""))
  res
}
