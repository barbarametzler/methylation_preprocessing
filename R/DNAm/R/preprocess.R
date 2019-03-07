preprocess = function(platform, input.folder, min.beads=3, detection=0.05, return.intensities = FALSE, return.snps.r=FALSE) {
  stopifnot(platform %in% c("hm450"))
  stopifnot(dir.exists(input.folder))
  stopifnot(min.beads > 0)
  stopifnot(detection > 0 & detection <= 1)

  # Read manifest
  probes=readRDS(paste0('illumina_manifests/', sprintf("%s_probes.rds", platform)))
  control.beads=readRDS(paste0('illumina_manifests/', sprintf("%s_controls.rds", platform)))

  # Separate probes by type
  inf1grn <- which(probes$type == "I-Grn")
  inf1red <- which(probes$type == "I-Red")
  inf2 <- which(probes$type == "II")

  idat.files <- list.idat(input.folder)
  idat.files$decoding <- as.POSIXct(NA)
  idat.files$scan <- as.POSIXct(NA)

  intensities <- array(NA, c(2, nrow(idat.files), nrow(probes)),
                       dimnames=list(c("A", "B"),
                                     idat.files$sample.id,
                                     rownames(probes)))

  controls <- array(NA, c(2, nrow(idat.files), nrow(control.beads)),
                    dimnames=list(c("grn", "red"),
                                  idat.files$sample.id,
                                  rownames(control.beads)))

  z <- qnorm(1 - detection)

  t0=Sys.time()
  for (i in 1:nrow(idat.files)) {
    print(paste0('Sample no. ', i, ' - time: ', Sys.time()))
    # Parse IDAT (function defined in illuminaio.R)
    data <- parse.idat(idat.files$grn[i], idat.files$red[i])
    stopifnot(data$barcode == idat.files$chip[i])
    stopifnot(data$label == idat.files$position[i])

    if (!is.null(data$run.info)) {
      idat.files$decoding[i] <- with(data$run.info, min(ts[type == "Decoding"]))
      idat.files$scan[i] <- with(data$run.info, min(ts[type == "Scan"]))
    }

    # Separate Grn/Red intensities into A (unmethylated) and B (methylated) intensities
    # depending on their type

    idx <- as.character(probes$address.a[inf1grn])
    intensities["A",idat.files$sample.id[i],inf1grn] <-
      ifelse(data$beads[idx,"grn.n"] >= min.beads, data$beads[idx,"grn.mean"], NA)

    idx <- as.character(probes$address.b[inf1grn])
    intensities["B",idat.files$sample.id[i],inf1grn] <-
      ifelse(data$beads[idx,"grn.n"] >= min.beads, data$beads[idx,"grn.mean"], NA)

    idx <- as.character(probes$address.a[inf1red])
    intensities["A",idat.files$sample.id[i],inf1red] <-
      ifelse(data$beads[idx,"red.n"] >= min.beads, data$beads[idx,"red.mean"], NA)

    idx <- as.character(probes$address.b[inf1red])
    intensities["B",idat.files$sample.id[i],inf1red] <-
      ifelse(data$beads[idx,"red.n"] >= min.beads, data$beads[idx,"red.mean"], NA)

    idx <- as.character(probes$address.a[inf2])
    intensities["A",idat.files$sample.id[i],inf2] <-
      ifelse(data$beads[idx,"red.n"] >= min.beads, data$beads[idx,"red.mean"], NA)
    intensities["B",idat.files$sample.id[i],inf2] <-
      ifelse(data$beads[idx,"grn.n"] >= min.beads, data$beads[idx,"grn.mean"], NA)

    idx <- as.character(rownames(control.beads))
    controls["grn",idat.files$sample.id[i],idx] <-
      ifelse(data$beads[idx,"grn.n"] > 0, data$beads[idx,"grn.mean"], NA)
    controls["red",idat.files$sample.id[i],idx] <-
      ifelse(data$beads[idx,"red.n"] > 0, data$beads[idx,"red.mean"], NA)


    # Define threshold of detection based on 600 negative control beads
    neg.beads <- which(control.beads$type == "NEGATIVE")
    neg.means <- rowMeans(controls[,idat.files$sample.id[i],neg.beads], na.rm=TRUE)
    neg.sds <- apply(controls[,idat.files$sample.id[i],neg.beads], 1, sd, na.rm=TRUE)
    threshold.inf1grn <- 2 * neg.means["grn"] + z * sqrt(2) * neg.sds["grn"]
    threshold.inf1red <- 2 * neg.means["red"] + z * sqrt(2) * neg.sds["red"]
    threshold.inf2 <- sum(neg.means) + z * sqrt(sum(neg.sds**2))

    # Background substraction
    I <- colSums(intensities[,idat.files$sample.id[i],])

    intensities["A",idat.files$sample.id[i],inf1grn] <- ifelse(
      intensities["A",idat.files$sample.id[i],inf1grn] > neg.means["grn"] &
        I[inf1grn] > threshold.inf1grn,
      intensities["A",idat.files$sample.id[i],inf1grn] - neg.means["grn"],
      NA
    )
    intensities["B",idat.files$sample.id[i],inf1grn] <- ifelse(
      intensities["B",idat.files$sample.id[i],inf1grn] > neg.means["grn"] &
        I[inf1grn] > threshold.inf1grn,
      intensities["B",idat.files$sample.id[i],inf1grn] - neg.means["grn"],
      NA
    )

    intensities["A",idat.files$sample.id[i],inf1red] <- ifelse(
      intensities["A",idat.files$sample.id[i],inf1red] > neg.means["red"] &
        I[inf1red] > threshold.inf1red,
      intensities["A",idat.files$sample.id[i],inf1red] - neg.means["red"],
      NA
    )
    intensities["B",idat.files$sample.id[i],inf1red] <- ifelse(
      intensities["B",idat.files$sample.id[i],inf1red] > neg.means["red"] &
        I[inf1red] > threshold.inf1red,
      intensities["B",idat.files$sample.id[i],inf1red] - neg.means["red"],
      NA
    )

    intensities["A",idat.files$sample.id[i],inf2] <- ifelse(
      intensities["A",idat.files$sample.id[i],inf2] > neg.means["red"] &
        I[inf2] > threshold.inf2,
      intensities["A",idat.files$sample.id[i],inf2] - neg.means["red"],
      NA
    )
    intensities["B",idat.files$sample.id[i],inf2] <- ifelse(
      intensities["B",idat.files$sample.id[i],inf2] > neg.means["grn"] &
        I[inf2] > threshold.inf2,
      intensities["B",idat.files$sample.id[i],inf2] - neg.means["grn"],
      NA
    )

    # Extract normalization probes for Grn and Red, and form the dye bias correction constant
    norm.grn.beads <- which(control.beads$type %in% c("NORM_C", "NORM_G"))
    norm.red.beads <- match(
      chartr("CG", "TA", control.beads$description[norm.grn.beads]),
      control.beads$description
    )
    norm.data <- cbind(
      grn=controls["grn",idat.files$sample.id[i],norm.grn.beads],
      red=controls["red",idat.files$sample.id[i],norm.red.beads]
    )
    corrections <- colMeans(
      rowMeans(norm.data, na.rm=TRUE) / norm.data,
      na.rm=TRUE
    )

    # Apply dye bias correction
    intensities["A",idat.files$sample.id[i],inf2] <-
      intensities["A",idat.files$sample.id[i],inf2] * corrections["red"]
    intensities["B",idat.files$sample.id[i],inf2] <-
      intensities["B",idat.files$sample.id[i],inf2] * corrections["grn"]
  }

  idx <- substr(rownames(probes), 1, 2) == "rs"

  # Create DNAm ratios as B (methylated) over total
  dnam <- intensities["B",,!idx] / colSums(intensities[,,!idx])

  # Extract intensities for SNP beads using theta/r format (as in GenomeStudio)
  #snps <- array(NA, c(2, nrow(idat.files), sum(idx)),
  #              dimnames=list(c("theta", "r"),
  #                            idat.files$sample.id,
  #                            rownames(probes)[idx]))
  #snps["theta",,] <-
  #  atan2(intensities["B",,idx], intensities["A",,idx]) / (pi / 2)
  #snps["r",,] <- sqrt(colSums(intensities[,,idx]**2))
  snps <- array(NA, c(nrow(idat.files), sum(idx)),
                dimnames=list(idat.files$sample.id,
                              rownames(probes)[idx]))
  # theta format
  snps <- atan2(intensities["B",,idx], intensities["A",,idx]) / (pi / 2)
  # r format
  if (return.snps.r) {
    snps_r <- sqrt(colSums(intensities[,,idx]**2))
  }


  # Extract all control probes data, and add summary statistics to samples table

  idx.bg <- match(c(
    "BS Conversion I-U1",
    "BS Conversion I-U2",
    "BS Conversion I-U3"
  ), control.beads$description)
  idx.signal <- match(
    chartr("U", "C", control.beads$description[idx.bg]),
    control.beads$description
  )
  idat.files$bc1.grn <- rowMeans(controls["grn",,idx.signal], na.rm=TRUE) /
    rowMeans(controls["grn",,idx.bg], na.rm=TRUE)

  idx.bg <- na.omit(match(c(
    "BS Conversion I-U4",
    "BS Conversion I-U5",
    "BS Conversion I-U6" # Removed in EPIC
  ), control.beads$description))
  idx.signal <- match(
    chartr("U", "C", control.beads$description[idx.bg]),
    control.beads$description
  )
  idat.files$bc1.red <- rowMeans(controls["red",,idx.signal], na.rm=TRUE) /
    rowMeans(controls["red",,idx.bg], na.rm=TRUE)

  idx <- which(control.beads$type == "BISULFITE CONVERSION II")
  idat.files$bc2 <- rowMeans(controls["red",,idx], na.rm=TRUE) /
    rowMeans(controls["grn",,idx], na.rm=TRUE)

  idat.files$ext.a <- controls["red",,control.beads$description == "Extension (A)"]
  idat.files$ext.c <- controls["grn",,control.beads$description == "Extension (C)"]
  idat.files$ext.g <- controls["grn",,control.beads$description == "Extension (G)"]
  idat.files$ext.t <- controls["red",,control.beads$description == "Extension (T)"]
  idat.files$hyb.low <- controls["grn",,control.beads$description == "Hyb (Low)"]
  idat.files$hyb.med <- controls["grn",,control.beads$description == "Hyb (Medium)"]
  idat.files$hyb.high <- controls["grn",,control.beads$description == "Hyb (High)"]
  idat.files$np.a <- controls["red",,control.beads$description == "NP (A)"]
  idat.files$np.c <- controls["grn",,control.beads$description == "NP (C)"]
  idat.files$np.g <- controls["grn",,control.beads$description == "NP (G)"]
  idat.files$np.t <- controls["red",,control.beads$description == "NP (T)"]

  idx.bg <- match(c(
    "GT Mismatch 1 (MM)",
    "GT Mismatch 2 (MM)",
    "GT Mismatch 3 (MM)"
  ), control.beads$description)
  idx.signal <- match(
    gsub("MM", "PM", control.beads$description[idx.bg], fixed=TRUE),
    control.beads$description
  )
  idat.files$spec1.grn <- rowMeans(controls["grn",,idx.signal], na.rm=TRUE) /
    rowMeans(controls["grn",,idx.bg], na.rm=TRUE)

  idx.bg <- match(c(
    "GT Mismatch 4 (MM)",
    "GT Mismatch 5 (MM)",
    "GT Mismatch 6 (MM)"
  ), control.beads$description)
  idx.signal <- match(
    gsub("MM", "PM", control.beads$description[idx.bg], fixed=TRUE),
    control.beads$description
  )
  idat.files$spec1.red <- rowMeans(controls["red",,idx.signal], na.rm=TRUE) /
    rowMeans(controls["red",,idx.bg], na.rm=TRUE)

  idx <- which(control.beads$type == "SPECIFICITY II")
  idat.files$spec2 <- rowMeans(controls["red",,idx], na.rm=TRUE) /
    rowMeans(controls["grn",,idx], na.rm=TRUE)

  idx.bg <- which(control.beads$description == "Biotin (Bkg)")
  idx.signal <- which(control.beads$description == "Biotin (High)")
  idat.files$st.grn <- controls["grn",,idx.signal] /
    controls["grn",,idx.bg]

  idx.bg <- which(control.beads$description == "DNP (Bkg)")
  idx.signal <- which(control.beads$description == "DNP (High)")
  idat.files$st.red <- controls["red",,idx.signal] /
    controls["red",,idx.bg]

  idx <- which(control.beads$type == "TARGET REMOVAL")
  idat.files$tr <- apply(controls["grn",,idx], 1, max, na.rm=TRUE)

  idat.files$missing <- rowMeans(is.na(dnam))

  idat.files$median.chrX <-
    apply(dnam[,rownames(subset(probes, chr == "X"))], 1, median, na.rm=TRUE)

  idat.files$missing.chrY <-
    rowMeans(is.na(dnam[,rownames(subset(probes, chr == "Y"))]))

  idat.files$grn <- NULL
  idat.files$red <- NULL

  # Save raw intensities (A/B as above) for all probes (including controls),
  # DNAm values (ratios), SNPs (theta/r format), and samples table
  #saveRDS(intensities,
  #        file.path(output, paste(name, "intensities.rds", sep="_")))
  #saveRDS(controls,
  #        file.path(output, paste(name, "controls.rds", sep="_")))
  #saveRDS(dnam,
  #        file.path(output, paste(name, "dnam.rds", sep="_")))
  #saveRDS(snps,
  #        file.path(output, paste(name, "snps.rds", sep="_")))
  #saveRDS(idat.files,
  #        file.path(output, paste(name, "samples.rds", sep="_")))

  t1=Sys.time()
  print(t1-t0)

  # Return list with sample table, SNPs theta values, CpG DNAm ratios and optionally: SNPs r values, intensities, controls
  rownames(idat.files) = idat.files$sample.id
  idat.files = idat.files[, 2:ncol(idat.files)]
  processed = list(
    samples = idat.files,
    cpgs = dnam,
    snps = snps
  )

  if (return.snps.r) {
    processed$snps_r = snps_r
  }

  if (return.intensities) {
    processed$intensities_A = intensities["A",,]
    processed$intensities_B = intensities["B",,]
    processed$controls_red = controls["red",,]
    processed$controls_grn = controls["grn",,]
  }

  processed
}
