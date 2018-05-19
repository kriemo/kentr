
#' Read in bam file as a trbl_interval
#' @param filename path to bam file
#' @export
bam_to_df <- function(filename = NULL,
                      region = ".",
                      tags = ""){
  if(is.null(filename)) filename <- system.file("extdata", "small_sorted.bam", package = "valr")

  filename <- path.expand(filename)

  res <- read_bam(filename, region, tags)

  res
}
