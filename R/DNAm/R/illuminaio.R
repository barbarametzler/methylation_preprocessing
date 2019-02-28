#tidyr::extract
#tidyr::spread

list.idat <- function(path) {
    file_paths <- list.files(path, "\\.idat$", full.names=TRUE, recursive=TRUE)

    files = tidyr::extract(
      data.frame(path=normalizePath(file_paths)),
      path,
      into = c("sample.id", "chip", "position", "channel"),
      regex = "(([0-9]+)_([RC0-9]+))_(Grn|Red)\\.idat$",
      remove=FALSE
    )
    files$path = as.character(files$path)
    files = tidyr::spread(files, channel, path)
    files = files[order(files$sample.id),]
    colnames(files) = tolower(colnames(files))

    files
}

parse.idat <- function(path.grn, path.red) {
    grn <- read.idat(path.grn)
    red <- read.idat(path.red)

    if (grn$format != red$format) {
        stop("Format does not match between Grn and Red files")
    }
    if (grn$barcode != red$barcode) {
        stop("Barcode does not match between Grn and Red files")
    }
    if (grn$label != red$label) {
        stop("Label does not match between Grn and Red files")
    }
    if (any(rownames(grn$beads) != rownames(red$beads))) {
        stop("Bead IDs do not match between Grn and Red files")
    }

    result <- list(
        format=grn$format,
        barcode=grn$barcode,
        label=grn$label,
        beads=cbind(grn$beads, red$beads)
    )
    colnames(result$beads) <- paste(rep(c("grn", "red"),
                                        c(ncol(grn$beads), ncol(red$beads))),
                                    colnames(result$beads),
                                    sep=".")

    if (!is.null(grn$run.info)) {
        result$run.info <- grn$run.info
    } else if (!is.null(red$run.info)) {
        result$run.info <- red$run.info
    }

    result
}

read.bpm <- function(input.file) {
    if (is.character(input.file)) {
        input.file <- file(input.file, "rb")
        on.exit(close(input.file))
    }
    if (!inherits(input.file, "connection")) {
        stop("`input.file` must be a character string or a connection")
    }
    if (!isOpen(input.file, "rb")) {
        open(input.file, "rb")
        on.exit(close(input.file))
    }

    # Read and check header
    header <- readChar(input.file, 3)
    if (header != "BPM") {
        stop("Invalid header")
    }

    # Read and check version
    version <- c(
        readBin(input.file, what="int", size=1, endian="little", signed=FALSE),
        readBin(input.file, what="int", size=4, endian="little", signed=TRUE)
    )
    if (!all.equal(version, c(1, 4))) {
        stop("Invalid version")
    }

    chip.type <- read.string(input.file)
    controls <- read.string(input.file)
    csv.lines <- readLines(input.file)

    csv.lines
}

read.idat <- function(input.file) {
    if (is.character(input.file)) {
        input.file <- file(input.file, "rb")
        on.exit(close(input.file))
    }
    if (!inherits(input.file, "connection")) {
        stop("`input.file` must be a character string or a connection")
    }
    if (!isOpen(input.file, "rb")) {
        open(input.file, "rb")
        on.exit(close(input.file))
    }

    result <- NULL

    # Read and check header
    header <- readChar(input.file, 4)
    if (header != "IDAT") {
        stop("Invalid header")
    }

    # Read and check version
    version <- readBin(input.file, what="int", n=2, size=4, endian="little")
    if (!all.equal(version, c(3, 0))) {
        stop("Invalid version")
    }
    result$version <- paste(version, collapse=".")

    # Read number of fields, field codes and offsets
    n.fields <- readBin(input.file, what="int", size=4, endian="little")
    fields <- matrix(NA, n.fields, 2, dimnames=list(NULL, c("code", "offset")))
    for (i in 1:n.fields) {
        fields[i,"code"] <- readBin(input.file, what="int", size=2,
                                    signed=FALSE, endian="little")
        fields[i,"offset"] <- readBin(input.file, what="int", size=8,
                                      endian="little")
    }

    # Read format
    i <- which(fields[,"code"] == 403)
    if (length(i) == 1) {
        seek(input.file, fields[i,"offset"])
        result$format <- read.string(input.file)
    }

    # Read barcode
    i <- which(fields[,"code"] == 402)
    if (length(i) == 1) {
        seek(input.file, fields[i,"offset"])
        result$barcode <- read.string(input.file)
    }

    # Read label
    i <- which(fields[,"code"] == 404)
    if (length(i) == 1) {
        seek(input.file, fields[i,"offset"])
        result$label <- read.string(input.file)
    }

    # Read run information
    i <- which(fields[,"code"] == 300)
    if (length(i) == 1) {
        seek(input.file, fields[i,"offset"])
        n.entries <- readBin(input.file, what="int", size=4, endian="little")
        if (n.entries > 0) {
            result$run.info <- matrix(NA, n.entries, 5,
                dimnames=list(NULL, c("ts", "type", "params", "software", "version")))
            for (i in 1:n.entries) for (j in 1:5) {
                result$run.info[i,j] <- read.string(input.file)
            }
            result$run.info <- as.data.frame(result$run.info, stringsAsFactors=FALSE)
            result$run.info$ts <- strptime(result$run.info$ts, "%m/%d/%Y %I:%M:%S %p")
        }
    }

    # Read number of beads
    i <- which(fields[,"code"] == 1000)
    if (length(i) != 1) {
        stop("Unable to determine number of beads")
    }
    seek(input.file, fields[i,"offset"])
    n.beads <- readBin(input.file, what="int", size=4, endian="little")
    stopifnot(n.beads > 0)

    # Allocate memory for beads
    result$beads <- matrix(NA, n.beads, 3,
                           dimnames=list(NULL, c("n", "mean", "sd")))

    # Read bead IDs
    i <- which(fields[,"code"] == 102)
    if (length(i) != 1) {
        stop("Unable to read bead IDs")
    }
    seek(input.file, fields[i,"offset"])
    rownames(result$beads) <- readBin(input.file, what="int", n=n.beads,
                                      size=4, endian="little")

    # Read bead counts
    i <- which(fields[,"code"] == 107)
    if (length(i) != 1) {
        stop("Unable to read bead counts")
    }
    seek(input.file, fields[i,"offset"])
    result$beads[,"n"] <- readBin(input.file, what="int", n=n.beads,
                                  size=1, signed=FALSE, endian="little")

    # Read intensity means
    i <- which(fields[,"code"] == 104)
    if (length(i) != 1) {
        stop("Unable to read intensity means")
    }
    seek(input.file, fields[i,"offset"])
    result$beads[,"mean"] <- readBin(input.file, what="int", n=n.beads,
                                     size=2, signed=FALSE, endian="little")

    # Read intensity SDs
    i <- which(fields[,"code"] == 103)
    if (length(i) != 1) {
        stop("Unable to read intensity SDs")
    }
    seek(input.file, fields[i,"offset"])
    result$beads[,"sd"] <- readBin(input.file, what="int", n=n.beads, size=2,
                                   signed=FALSE, endian="little")

    result
}

read.string <- function(con) {
    m <- readBin(con, what="int", size=1, signed=FALSE, endian="little")
    n <- m %% 128
    shift <- 0L
    while (m %/% 128 == 1) {
        m <- readBin(con, what="int", size=1, signed=FALSE, endian="little")
        shift <- shift + 7L
        k <- (m %% 128) * 2**shift
        n <- n + k
    }
    readChar(con, n)
}

