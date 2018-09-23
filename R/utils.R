
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

#' Download a file if it doesn't exist
#' @param url download url
#' @param filename output file name
#' @param gunzip uncompress file using gzip Default = FALSE, if suffix contains .gz this will
#' be removed.
#' @export
download_file <- function(url, filename, gunzip = FALSE){
  if(!file.exists(filename)){
    utils::download.file(url, filename)
  }

  if(gunzip){
    if(endsWith(filename, ".gz")) {
      R.utils::gunzip(filename)
    } else {
      stop(paste0(filename, " not a gzipped file"))
    }
  }
}

#' Read multiple tsv files and join by common cols
#' @param .dir directory path
#' @param .pattern wildcard for file matching
#' @param .cols joins to join multiple files
#' @param ... additional arguments to pass to [readr::read_tsv()]
#' @export
n_inner_join <- function(.dir, .pattern, .cols, ...){
  # read in and aggregate data into flat format
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

