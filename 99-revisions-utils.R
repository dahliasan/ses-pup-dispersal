## internal rotate to match old behaviour
## https://r-forge.r-project.org/scm/viewvc.php/pkg/raster/R/rotate.R?root=raster&r1=2782&r2=2981
#' importFrom raster merge
.rotate <- function(x, ...) {
    e <- raster::extent(x)
    xrange <- e@xmax - e@xmin
    inverse <- FALSE
    if (xrange < 350 | xrange > 370 | e@xmin < -10 | e@xmax > 370) {
        if (xrange < 350 | xrange > 370 | e@xmin < -190 | e@xmax > 190) {
            warning("this does not look like an appropriate object for this function")
        } else {
            inverse <- TRUE
        }
    }
    hx <- e@xmin + xrange / 2
    r1 <- crop(x, extent(e@xmin, hx, e@ymin, e@ymax))
    r2 <- crop(x, extent(hx, e@xmax, e@ymin, e@ymax))
    if (inverse) {
        r1@extent@xmin <- r2@extent@xmax
        r1@extent@xmax <- r1@extent@xmin + 0.5 * xrange
    } else {
        r2@extent@xmin <- r2@extent@xmin - xrange
        r2@extent@xmax <- r2@extent@xmax - xrange
    }
    ln <- names(x)
    out <- merge(r1, r2, overlap = FALSE, ...)
    names(out) <- names(x)
    out@z <- x@z

    # suggested by Mike Sumner:
    p <- projection(out)
    if (length(grep("\\+over", p)) > 0) {
        projection(out) <- gsub("[[:space:]]\\+over", "", p)
    }

    return(out)
}
