#' read a whole genome into memory as a dataframe
#'
#' @param genome_file path to fasta file that is indexed with samtools faidx
#'
#' @examples
#' fa_path <- system.file("extdata", "test.fasta", package = "kentr")
#' read_genome_seq(fa_path)
#'
#' @export
read_genome_seq <- function(genome_file) {
  gnome_df <- suppressWarnings(readr::read_tsv(paste0(genome_file,
                                                      ".fai"),
                                               col_names = c("chrom", "end"),
                                               col_types = c('ci')))

  gnome_df <- gnome_df %>%
    mutate(start = 0) %>%
    select(chrom, start, end)

  res <- get_sequences(gnome_df, genome_file)

  res <- mutate(res,
                seq = stringr::str_to_upper(seq))

  res
}

#' Find positions of query seq matches in seq column
#'
#' @param df data frame containing a header and sequence column
#' @param query_seq sequence to identify in df$seq column, can contain regex
#'
#' @return A data.frame with chrom, start, end, and header columns.
#' Start and end indicate start
#' point and end point of query_seq match.
#'
#'
#' @examples
#' fa_path <- system.file("extdata", "test.fasta", package = "kentr")
#' query_seq <- "AATAAA[GTC]"
#'
#' read_genome_seq(fa_path) %>%
#'   find_motifs(query_seq)
#'
#' @export
find_motifs <- function(sequence_df, query_seq) {

  if(!all(c("seq", "header") %in% colnames(sequence_df))){
    stop("columns named seq and header necessary in dataframe",
         call. = FALSE)
  }

  res <- stringr::str_locate_all(sequence_df[["seq"]],
                                 query_seq)
  names(res) <- sequence_df[["header"]]

  res <- purrr::imap_dfr(res,
                         function(indexes, seq_header) {
                           tibble::as_tibble(indexes) %>%
                             mutate(
                               start = start - 1,
                               header = seq_header,
                               chrom = stringr::str_remove(header, ":.+$")
                             ) %>%
                             select(chrom, start, end, header)
                         })

  res
}

