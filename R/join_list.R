#' Title Combind a list of named vector into a data.frame by name
#'
#' @param list a list of named vector, which each element has common name.
#'
#' @return a data frame
#'
#'
#'
join_list <- function(list){
  one_row_dataframe_list <- lapply(list,function(each_named_vector){
    each_named_vector <- as.data.frame(each_named_vector)
    each_named_vector <- tibble::rownames_to_column(each_named_vector,var = "term")
    return(each_named_vector)
  })
  #perform t() to
  result <- Reduce(function(x,y) merge(x,y,by="id",all=TRUE),one_row_dataframe_list, accumulate = FALSE)
  return(t(result))
}
