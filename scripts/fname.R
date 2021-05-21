#' Get a filename from the config file
#'
#' Convenience function for getting a filename from the config file.
#'
#' @param alias The short alias of the filename you want.
#' @param cell_line When specified, the cell line specific file will be returned.
#'
#' @return filename The requested filename.
#'
#' @examples
#' fname('annotations')
#' fname('profiles', cell_line = 'K652')
fname <- function(alias, cell_line = '') {
    config <- read_yaml('workflow/config.yaml')
    if (cell_line == '')
        paste0(config$data_dir, '/', config$filenames[[alias]])
    else
        paste0(config$data_dir, '/', cell_line, '/', config$filenames[[cell_line]][[alias]])
}
