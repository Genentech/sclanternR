
.log_msg <- function(level, verbose, fmt, ...) {
    if (verbose >= level) {
        message(sprintf(paste0("[%s] ", fmt), format(Sys.time(), "%H:%M:%S"), ...))
    }
}
