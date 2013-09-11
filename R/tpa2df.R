#' Function to convert TIMESAT .tpa binary format file to an R dataframe.
#'
#' @export
#' @param tpa_file_name A string giving the location of a .tpa file output by 
#' TIMESAT
#' @return A data.frame containing 14 columns: row, col, season, start, end, 
#' length, base_value, peak_time, peak_value, amp, left_deriv, right_deriv, 
#' large_integ, and small_integ
#' @examples
#' # TODO: Need to add examples here, and need to include a sample TIMESAT tpa 
#' # file in the package data.
tpa2df <- function(tpa_file_name) {
    require(base) # Needed for file.info

    if (missing(tpa_file_name) || !grepl('[.]tpa$', tolower(tpa_file_name))) {
        stop('must specify a .tpa file')
    }

    # The number of seasonal indicators output by TIMESAT.
    SEASONAL_INDICATORS <- 11

    # Number of elements in the tpa file line header (which are normally: row, 
    # column, number of seasons).
    LINE_HEADER_SIZE <- 3

    tpa_file_obj <- file(tpa_file_name, "rb")
    raw_vector <- readBin(tpa_file_obj, n=file.info(tpa_file_name)$size, raw())
    close(tpa_file_obj)

    # This function is used to track the offset within the binary vector as readBin 
    # does not track position except for file objects
    offset <- 1
    raw_vec_length <- length(raw_vector)
    offset_readBin <- function(raw_vec, what, n=n, size=size, increment_offset=TRUE, ...) {
        bin_data <- readBin(raw_vec[offset:(n * size + offset)], what, n, size, ...)
        # Use a global variable to track the offset
        if (increment_offset) {assign("offset", offset + (size*n), inherits=TRUE)}
        return(bin_data)
    }

    # File header format is: nyears nptperyear rowstart rowstop colstart colstop
    file_header <- offset_readBin(raw_vector, integer(), n=6, size=4)
    num_years <- file_header[1]
    rowstart <- file_header[3]
    rowstop <- file_header[4]
    colstart <- file_header[5]
    colstop <- file_header[6]

    num_pixels <- (colstop - colstart) * (rowstop - rowstart)

    # Setup tpa_data matrix to hold the data. Need to read first line of file 
    # to get number of seasons in order to set data matrix dimensions.
    first_line_header <- offset_readBin(raw_vector, integer(), n=3, size=4, increment_offset=FALSE)
    num_seasons <- first_line_header[3]
    tpa_data <- matrix(nrow=num_pixels*num_seasons,
                       ncol=(LINE_HEADER_SIZE + SEASONAL_INDICATORS))

    # Read the data and save it in the tpa_data matrix
    for (pixelnum in 1:num_pixels) {
        line_header <- offset_readBin(raw_vector, integer(), n=3, size=4)
        # Line header format is: rownum colnum numseasons
        if (line_header[3] != num_seasons) {
            stop(paste('pixel', pixelnum, 'has', line_header[3],
                       'seasons, but first line of .tpa file had', num_seasons, 
                       'seasons'))
        }
        for (seasonnum in 1:num_seasons) {
            line_data <- offset_readBin(raw_vector, numeric(), 
                                        n=SEASONAL_INDICATORS, size=4)
            tpa_data_row <- (pixelnum-1)*num_seasons + seasonnum
            tpa_data[tpa_data_row, ] <- c(line_header[1], line_header[2],
                                          seasonnum, line_data)
        }
    }

    tpa_data <- data.frame(tpa_data)
    names(tpa_data) <- c("row", "col", "season", "start", "end", "length",
                         "base_value", "peak_time", "peak_value", "amp", "left_deriv",
                         "right_deriv", "large_integ", "small_integ")
    return(tpa_data)
}
