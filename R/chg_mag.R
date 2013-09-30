#' Change Magnitude Image for CVAPS
#'
#' This code calculate the change magnitude image for the Change Vector 
#' Analysis in Posterior Probability Space (CVAPS) method of Chen et al. 2011.  
#' Use the change magnitude image in conjunction with the change direction 
#' image from \code{chg_dir}, and \code{DFPS} to use the Double Window Flexible 
#' Pace Search method (Chen et al. 2003) to determine the threshold to use to 
#' map areas of change and no-change.
#'
#' @export
#' @param x time 0 posterior probability \code{Raster*}
#' @param y time 1 posterior probability \code{Raster*}
#' @param filename (optional) filename for output change magnitude 
#' \code{RasterLayer}
#' @param \code{Raster*} object with change magnitude image
#' @references Chen, J., P. Gong, C. He, R. Pu, and P. Shi. 2003.
#' Land-use/land-cover change detection using improved change-vector analysis.
#' Photogrammetric Engineering and Remote Sensing 69:369-380.
#' 
#' Chen, J., X. Chen, X. Cui, and J. Chen. 2011. Change vector analysis in 
#' posterior probability space: a new method for land cover change detection.  
#' IEEE Geoscience and Remote Sensing Letters 8:317-321.
chg_mag <- function(x, y, filename="", ...) {
    out <- raster(x)
    if (proj4string(x) != proj4string(y)) {
        stop('Error: t0 and t1 coordinate systems do not match')
    }
    if (extent(x) != extent(y)) {
        stop('Error: t0 and t1 extents do not match')
    }
    if (nlayers(x) != nlayers(y)) {
        stop('Error: t0 and t1 probability maps have differing number of classes')
    }

    n_classes <- nlayers(x)

    cl <- getCluster()
    on.exit(returnCluster())

    nodes <- length(cl)

    # Process over blocks to conserve memory and use multiple cores.
    bs <- blockSize(x, minblocks=nodes*4)
    pb <- pbCreate(bs$n)

    calc_chg_mag <- function(i) {
        x_block <- getValues(x, row=bs$row[block_num], 
                             nrows=bs$nrows[block_num])
        y_block <- getValues(y, row=bs$row[block_num], 
                             nrows=bs$nrows[block_num])
        # Calculate change magnitude (eqn 3 in Chen 2011)
        dP <- y_block - x_block
        chgmag <- sqrt(rowSums((y_block - x_block)^2))
        return(chgmag)
    }

    # Get all nodes going.
    for (i in 1:nodes) {
        sendCall(cl[[i]], calc_chg_mag, i, tag=i)
    }

    # Save to temp file if cannot hold result in memory.
    filename <- trim(filename)
    if (!canProcessInMemory(out) & filename == "") {
        filename <- rasterTmpFile()
    }
    if (filename != "") {
        out <- writeStart(out, filename=filename, ... )
    } else {
        vv <- matrix(ncol=nrow(out), nrow=ncol(out))
    }

    pb <- pbCreate(bs$n)
    for (i in 1:bs$n) {
        # receive results from a node
        d <- recvOneData(cl)

        # error?
        if (! d$value$success) {
            stop('cluster error')
        }

        # which block is this?
        b <- d$value$tag
        cat('received block: ', b, '\n')
        flush.console()

        if (filename != "") {
            out <- writeValues(out, d$value$value, bs$row[b])
        } else {
            cols <- bs$row[b]:(bs$row[b] + bs$nrows[b] - 1)
            vv[,cols] <- matrix(d$value$value, nrow=out@ncols)
        }

        # need to send more data?
        ni <- nodes + i
        if (ni <= bs$n) {
            sendCall(cl[[d$node]], calc_chg_mag, ni, tag=ni)
        }
        pbStep(pb)
    }

    if (filename != "") {
        out <- writeStop(out)
    } else {
        out <- setValues(out, as.vector(vv))
    }
    pbClose(pb)

    return(out)
}
