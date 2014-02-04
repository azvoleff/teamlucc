#' @importFrom bitops cksum
verify_download <- function(espa_url, local_path) {
    cksum_file <- tempfile()
    ret_code <- download.file(gsub('\\.tar\\.gz$', '.cksum', espa_url), 
                              cksum_file, mode="w", quiet=TRUE)
    if (ret_code != 0) {
        message(paste('Warning: problem downloading cksum for', local_path))
        return(1)
    } else {
        # TODO: Check return code and handle accordingly
        # The first element is the cksum, second is the expected file size in 
        # bytes, and the third is the filename
        espa_checksum <- scan(cksum_file, what=c('integer', 'integer', 
                                                 'character'), quiet=TRUE)
        unlink(cksum_file)
        local_size <- file.info(local_path)$size
        # TODO: Figure out how to compute a checksum in R that matches the checksum 
        # output ESPA gives. It appears the ESPA checksum is a CRC from 'cksum' 
        # command run on Linux. This is not a CRC-32 checksum, so the R digest 
        # package won't work for computing it. bitops has a function, 'cksum' that 
        # might work.
        #local_crc <- strtoi(digest(, algo="crc32", file=TRUE), base=16L)

        # f = file(local_path,"rb")
        # local_crc <- cksum(rawToChar(readBin(f, raw(), n=local_size)))
        # close(f)

        # if (espa_checksum[1] != local_crc) {
        #     return(2)
        # } else if (espa_checksum[2] != local_size) {
        if (espa_checksum[2] != local_size) {
            return(3)
        } else {
            return(0)
        }
    }
}

download_ESPA_file <- function(espa_url, output_path) {
    ret_code <- download.file(espa_url, output_path, mode="wb")
    if (ret_code != 0) {
        message(paste('Warning: problem downloading', output_path))
        return(1)
    } else if (verify_download(espa_url, output_path) != 0) {
        message(paste("Warning: checksum mismatch on", output_path))
        return(2)
    } else {
        return(0)
    }
}

#' Download a completed ESPA order
#'
#' Function to download a set of Landsat images from ESPA given a valid order 
#' ID and the email that placed the ESPA order.
#' 
#' @export
#' @importFrom stringr str_extract
#' @param email address used to place the order
#' @param order_ID the ESPA order ID
#' @param output_folder the folder to save output data in
#' @return used for the side effect of downloading Landsat scenes
espa_download <- function(email, order_ID, output_folder) {
    email_re <- '^[A-Z0-9._%+-]+@[A-Z0-9.-]+\\.[A-Z]{2,4}$'
    if (!grepl(email_re, email, ignore.case=TRUE)) {
        stop(paste(email, 'does not appear to be a valid email address'))
    }
    # Below is for old format ESPA order IDs (through end of 2013)
    #if (!grepl('^[0-9]{13}$', order_ID)) {
    if (!grepl('^[0-9]{6}-[0-9]{6}$', order_ID)) {
        stop(paste(order_ID, 'does not appear to be a valid ESPA order ID'))
    }
    if (!file_test('-d', output_folder)) {
        stop(paste(output_folder, 'does not appear to be a valid directory'))
    }

    # Parse ESPA page for download links
    email_noat <- gsub('@', '%40', email)
    espa_page <- scan(paste0("http://espa.cr.usgs.gov/status/", email_noat, 
                             "-", order_ID), what='character', quiet=TRUE)
    url_re <- paste0('http://espa\\.cr\\.usgs\\.gov/orders/', email, '-', 
                     order_ID, '/L[ET][0-9]{14}-SC[0-9]{14}\\.tar\\.gz')
    espa_urls <- espa_page[grepl(url_re, espa_page)]
    espa_urls <- str_extract(espa_urls, url_re)

    successes <- 0
    failures <- 0
    skips <- 0
    message(paste('Found', length(espa_urls), 'ESPA downloads.'))
    for (n in 1:length(espa_urls)) {
        espa_url <- espa_urls[n]
        img_file <- basename(espa_url)
        output_path <- file.path(output_folder, img_file)
        if (file.exists(output_path)) {
            if (verify_download(espa_url, output_path)) {
                message(paste(img_file, 'exists but has bad checksum - re-downloading file'))
            } else {
                message(paste(img_file, 'exists and has good checksum - skipping download'))
                skips <- skips + 1
                next
            }
        }
        if (download_ESPA_file(espa_url, output_path) == 0) {
            successes <- successes + 1
        } else {
            failures <- failures + 1
        }

    }
    message(paste(successes, "file(s) succeeded,", skips, "file(s) skipped,", 
                failures, "file(s) failed."))
}
