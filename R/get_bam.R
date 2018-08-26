
#' Read in bam file as a trbl_interval
#' @param filename path to bam file
#' @param region samtools region query
#' @param tags bam tags to return, supplied with type suffix, i.e. "XS:A" for character,
#' "AS:i" for integer, and "CB:Z" for string
#' @export
bam_to_df <- function(filename = NULL,
                      region = ".",
                      tags = NULL){
  if(is.null(filename)) filename <- system.file("extdata", "small_sorted.bam", package = "kentr")

  filename <- path.expand(filename)


  if(!is.null(tags)){
    #parse tag type
    if(any(!str_detect(tags, ":[ZiA]$"))){
      stop("missing type value in supplied tags, need :Z, :i, or :A suffix on tags")
    }
    sp_tags <- stringr::str_split(tags, ":", simplify = T)
    tag_types <- sp_tags[, 2]
    tag_ids <- sp_tags[, 1]

  } else{
    tag_types = ""
    tag_ids = ""
  }
  res <- read_bam(filename, region, tag_ids, tag_types)

  res
}
