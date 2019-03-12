#mvoutlier::pcout
#quadprog::solve.QP

remove_unreliable_samples = function(dnam, thresholdNA = 0.1) {
  # Remove samples with too many missing values
  dnam$samples = dnam$samples[dnam$samples$missing < thresholdNA, ]

  # Remove potential multivariate outliers according to BC controls
  bc1_grn_lower_threshold = mean(dnam$samples$bc1.grn) - 2 * sd(dnam$samples$bc1.grn)
  bc1_grn_upper_threshold = mean(dnam$samples$bc1.grn) + 2 * sd(dnam$samples$bc1.grn)
  bc1_red_lower_threshold = mean(dnam$samples$bc1.red) - 2 * sd(dnam$samples$bc1.red)
  bc1_red_upper_threshold = mean(dnam$samples$bc1.red) + 2 * sd(dnam$samples$bc1.red)
  bc2_lower_threshold = mean(dnam$samples$bc2) - 2 * sd(dnam$samples$bc2)
  bc2_upper_threshold = mean(dnam$samples$bc2) + 2 * sd(dnam$samples$bc2)

  keep = (dnam$samples$bc1.grn >  bc1_grn_lower_threshold) &
         (dnam$samples$bc1.grn <= bc1_grn_upper_threshold) &
         (dnam$samples$bc1.red >  bc1_red_lower_threshold) &
         (dnam$samples$bc1.red <= bc1_red_upper_threshold) &
         (dnam$samples$bc2 >  bc2_lower_threshold) &
         (dnam$samples$bc2 <= bc2_upper_threshold)

  #tech.vars = scale(dnam$samples[,c("bc1.grn", "bc1.red", "bc2")])
  #outliers <- mvoutlier::pcout(tech.vars)
  #dnam$samples = dnam$samples[outliers$wfinal01 == 1, ]
  dnam$samples = dnam$samples[keep, ]

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
    plot(dnam$samples$median.chrX, dnam$samples$missing.chrY, col=addNA(sex), xlim=c(0,1), ylim=c(0,1),
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

  calls[dnam$snps <= non.carrier.threshold] = 0 # non-carrier
  calls[dnam$snps > non.carrier.threshold & dnam$snps < homozygous.threshold] = 1 # heterozygous
  calls[dnam$snps >= homozygous.threshold] = 2 # homozygous

  if (plot) {
    title = paste(dim(dnam$snps)[1], "samples, ", dim(dnam$snps)[2], "SNPs")
    boxplot(as.numeric(as.matrix(dnam$snps)) ~ as.numeric(as.matrix(call_snps(dnam))), xlab="Carrier status", ylab="Theta intensities", main = title)
    abline(h=non.carrier.threshold)
    abline(h=homozygous.threshold)
  }

  calls
}

snp_distance = function(dnam, plot=FALSE) {
  dist = dist(call_snps(dnam), method="manhattan")

  if (plot) {
    hist(dist, main="Genotype distance among 65 SNPs", xlab="Absolute number of different alleles", freq=F)
  }

  dist
}

identify_replicates = function(dnam, snp.distance.error.margin = 0) {
  inferred.sex = infer_sex(dnam)
  snp.distance = snp_distance(dnam)

  # Identify possible replicate pairs: those with (near) zero SNP distance
  replicate_ids_in_lower_tri = which(snp.distance < snp.distance.error.margin)

  if (!length(replicate_ids_in_lower_tri)) return(0)

  # Initalize empty dataframe (potentially too large)
  replicates = data.frame(a = character(length(replicate_ids_in_lower_tri)),
                          b = character(length(replicate_ids_in_lower_tri)),
                          stringsAsFactors = F)

  # Iterate over all possible pairs
  for (i in 1:length(replicate_ids_in_lower_tri)) {
    # Identify pair from position in distance matrix
    # https://stackoverflow.com/questions/30376553/get-labels-of-distance-matrix-cell
    pair = labels(snp.distance)[which(lower.tri(snp.distance),arr.ind=TRUE)[replicate_ids_in_lower_tri[i],]]

    # Make sure sex is the same
    if (!is.na(inferred.sex[pair[1]]) & !is.na(inferred.sex[pair[2]])) {
      if (inferred.sex[pair[1]] == inferred.sex[pair[2]]) {
        replicates[i,] = pair
      }
    }
  }

  # Shrink dataframe to final size
  replicates = replicates[replicates$a != "",]
  rownames(replicates) = 1:nrow(replicates)

  replicates
}

estimate_leukocytes = function(dnam) {
  coefs = readRDS('illumina_manifests/leukocyte_coefs.rds')

  probe.ids <- intersect(colnames(dnam$cpgs), rownames(coefs))

  betas <- t(dnam$cpgs[,probe.ids])
  coefs <- coefs[probe.ids,]

  wbc.predictions <- matrix(NA, ncol(betas), ncol(coefs),
                            dimnames=list(colnames(betas), colnames(coefs)))
  A <- diag(ncol(coefs))
  b <- rep(0, ncol(coefs))
  for (i in 1:nrow(wbc.predictions)) {
    idx <- which(!is.na(betas[,i]))
    D <- crossprod(coefs[idx,])
    d <- crossprod(coefs[idx,], betas[idx,i])
    wbc.predictions[i,] <- quadprog::solve.QP(D, d, A, b)$solution
  }

  wbc.predictions
}
