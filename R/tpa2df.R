#' Function to convert TIMESAT .tpa binary format file to an R dataframe.
#'
#' @export
#' @param x A string giving the location of a .tpa file output by 
#' TIMESAT
#' @param max_num_seasons the maximum number of seasons for any of the pixels 
#' in the file
#' @return A data.frame containing 14 columns: row, col, season, start, end, 
#' length, base_value, peak_time, peak_value, amp, left_deriv, right_deriv, 
#' large_integ, and small_integ
#' @examples
#' # TODO: Need to add examples here, and need to include a sample TIMESAT tpa 
#' # file in the package data.
tpa2df <- function(x, max_num_seasons) {
    if (missing(x) || !grepl('[.]tpa$', tolower(x))) {
        stop('must specify a .tpa file')
    }
    if (missing(max_num_seasons) || max_num_seasons < 1) {
        stop('must specify maximum number of seasons represented in tpa file')
    }

    # The number of seasonal indicators output by TIMESAT.
    NUM_SEASONAL_INDICATORS <- 11

    # Number of elements in the tpa file line header (which are normally: row, 
    # column, number of seasons).
    LINE_HEADER_SIZE <- 3

    tpa_file_obj <- file(x, "rb")
    raw_vector <- readBin(tpa_file_obj, n=file.info(x)$size, raw())
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

    num_rows <- rowstop - rowstart + 1
    num_cols <- colstop - colstart + 1

    num_pixels <- num_cols * num_rows

    tpa_data <- matrix(nrow=num_pixels*max_num_seasons,
                       ncol=(LINE_HEADER_SIZE + NUM_SEASONAL_INDICATORS))

    # Read the data and save it in the tpa_data matrix
    for (pixelnum in 1:num_pixels) {
        line_header <- offset_readBin(raw_vector, integer(), n=3, size=4)
        # Line header format is: rownum colnum numseasons
        num_seasons <- line_header[3]
        if (num_seasons > max_num_seasons) {
            stop(paste('pixel', pixelnum, 'has', num_seasons,
                       'seasons, but max_num_seasons was set to ', max_num_seasons, 
                       'seasons'))
        }
        if (num_seasons == 0) {
            # No seasons identified for this pixel - skip
            next
        }
        for (seasonnum in 1:max_num_seasons) {
            # Seasons were found for this pixel - read them
            tpa_data_row <- (pixelnum - 1)*num_seasons + seasonnum
            line_data <- offset_readBin(raw_vector, numeric(), 
                                        n=NUM_SEASONAL_INDICATORS, size=4)
            tpa_data[tpa_data_row, ] <- c(line_header[1], line_header[2], 
                                          seasonnum, line_data)
        }
    }

    tpa_data <- data.frame(tpa_data)
    tpa_data <- tpa_data[!(rowSums(is.na(tpa_data)) == ncol(tpa_data)), ]
    names(tpa_data) <- c("row", "col", "season", "start", "end", "length",
                         "base_value", "peak_time", "peak_value", "amp", "left_deriv",
                         "right_deriv", "large_integ", "small_integ")
    return(tpa_data)
}
