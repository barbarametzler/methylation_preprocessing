#mvoutlier::pcout

remove_unreliable_samples = function(dnam, thresholdNA = 0.1) {
  # Remove samples with too many missing values
  dnam$samples = dnam$samples[dnam$samples$missing < thresholdNA, ]

  # Remove potential multivariate outliers according to BC controls
  tech.vars = scale(dnam$samples[,c("bc1.grn", "bc1.red", "bc2")])
  outliers <- mvoutlier::pcout(tech.vars)
  dnam$samples = dnam$samples[outliers$wfinal01 == 1, ]

  # Also remove corresponding samples from other list-entries in dnam (ratios, snps, possibly controls and intensities)
  dnam$cpgs = dnam$cpgs[rownames(dnam$samples), ]
  dnam$snps = dnam$snps[rownames(dnam$samples), ]
  if (!is.null(dnam$intensities_A)) dnam$intensities_A = dnam$intensities_A[rownames(dnam$samples), ]
  if (!is.null(dnam$intensities_B)) dnam$intensities_B = dnam$intensities_B[rownames(dnam$samples), ]
  if (!is.null(dnam$controls_red)) dnam$controls_red = dnam$controls_red[rownames(dnam$samples), ]
  if (!is.null(dnam$controls_grn)) dnam$controls_grn = dnam$controls_grn[rownames(dnam$samples), ]

  dnam
}


infer_sex = function(dnam, plot=F, threshold.chrX = 0.4, threshold.chrY = 0.4) {
  sex = factor(rep(NA, times=nrow(dnam$samples)), levels=c("f", "m"))
  names(sex) = rownames(dnam$sample)

  # Female
  sex[rownames(dnam$samples)[dnam$samples$median.chrX > threshold.chrX &
                             dnam$samples$missing.chrY > threshold.chrY]] = "f"

  # Male
  sex[rownames(dnam$samples)[dnam$samples$median.chrX <= threshold.chrX &
                               dnam$samples$missing.chrY <= threshold.chrY]] = "m"

  if (plot) {
    plot(dnam$samples$median.chrX, dnam$samples$missing.chrY, col=sex, xlim=c(0,1), ylim=c(0,1),
      xlab = "Median methylation level Chr. X",
      ylab = "Proportion missing Chr. Y")
    rect(threshold.chrX, threshold.chrY, 1, 1, border = "red")
    rect(0.0, 0.0, threshold.chrX, threshold.chrY, border = "blue")
  }

  sex
}


call_snps = function(dnam, non.carrier.threshold = 0.2, homozygous.threshold = 0.8, plot = FALSE) {
  calls = array(NA, dim=c(nrow(dnam$samples), ncol(dnam$snps)),
                dimnames=list(rownames(dnam$samples), colnames(dnam$snps)))

  calls[dnam$snps < non.carrier.threshold] = 0 # non-carrier
  calls[dnam$snps > non.carrier.threshold & dnam$snps < homozygous.threshold] = 1 # heterozygous
  calls[dnam$snps > homozygous.threshold] = 2 # homozygous

  if (plot) {
    title = paste(dim(eira$snps)[1], "samples, ", dim(eira$snps)[2], "SNPs")
    boxplot(as.numeric(as.matrix(eira$snps)) ~ as.numeric(as.matrix(call_snps(eira))), xlab="Carrier status", ylab="Theta intensities", main = title)
    abline(h=non.carrier.threshold)
    abline(h=homozygous.threshold)
  }

  calls
}


identify_replicates = function(dnam) {
  samples = dnam$samples
  samples$inferred.sex = infer_sex(dnam)
  calls = call_snps(dnam)

  # Merge SNPs
  samples <- cbind(samples, calls)

  # Compute similarity matrix
  K <- matrix(0, nrow(samples), nrow(samples),
              dimnames=list(rownames(samples), rownames(samples)))

  tmp <- t(data.matrix(samples[,c("inferred.sex", colnames(calls))]))

  for (i in 1:nrow(K)) {
    K[i,] <- as.integer(apply((tmp[,i] - tmp) == 0, 2, mean, na.rm=TRUE))
  }

  # Identify duplicates from SNPs
  dups <- as.data.frame(which(K == 1, arr.ind=TRUE))
  dups$sample1 <- rownames(samples)[dups$row]
  dups$sample2 <- rownames(samples)[dups$col]
  dups$id1 <- rownames(samples)[dups$row]
  dups$id2 <- rownames(samples)[dups$col]
  dups$source <- "snps"
  dups$row = dups$col <- NULL

  dups <- subset(dups,
                 sample1 != sample2
  )

  dups
}
