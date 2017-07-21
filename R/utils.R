#' @export

#' @export
rev_comp <- function(x, rev_complement = T){
  bps <- stringr::str_split(x, "", simplify = T)
  res <- purrr::map_chr(bps, bp_comp)
  if (rev_complement == T) res <- rev(res)
  res
}

bp_comp <- function(nt){
  #return complementary basepair
  comp <- c("A", "G", "T", "C", "N")
  names(comp) <- c("T", "C", "A", "G", "N")
  comp[nt] %>% unname()
}

#' @export
n_inner_join <- function(.dir, .pattern, .cols, ...){
  # read in and aggregate data into matrix-like format
  files <- dir(.dir,
               full.names = T,
               pattern = .pattern)
  dat <- suppressMessages(purrr::map(files,
                                     ~readr::read_tsv(.x,
                                                      progress = F,
                                                      ...
                                                      )))
  dat <- Reduce(function(x, y) dplyr::inner_join(x, y, by = .cols), dat)
  colnames(dat) <- c(.cols, basename(files))
  dat
}

#' @export
write_gztsv <- function(df, name, ...){
  # write output as gzipped, supply name without .gz
  if(stringr::str_detect(name, ".gz$")){
    uncompressed_name <- stringr::str_replace(name, ".gz$", "")
  } else {
    uncompressed_name <- name
  }
  readr::write_tsv(df, uncompressed_name, ...)
  system(paste0("gzip -f ", uncompressed_name))
}
