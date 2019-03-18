#quadprog::solve.QP

remove_unreliable_samples_probes = function(dnam, threshold_NA_sample = 0.1, threshold_NA_cpg = 0.2, bc_threshold = 1, plot=T) {
  # Define samples with failed bisulfite conversion (Illumina recommends threshold of 1 (Heiss 2019, https://doi.org/10.1186/s13148-019-0615-3)
  bc1_grn_failed = dnam$samples$bc1.grn < bc_threshold
  bc1_red_failed = dnam$samples$bc1.red < bc_threshold
  bc2_failed = dnam$samples$bc2 < bc_threshold

  # Define samples with too many missing values
  too_many_NAs_per_sample = dnam$samples$missing > threshold_NA_sample
  keep_samples = !bc1_grn_failed & !bc1_red_failed & !bc2_failed & !too_many_NAs_per_sample

  # Define CpG sites with too many missing values
  NA_proportion_columns = colMeans(is.na(dnam$cpgs[keep_samples,]), na.rm=T)
  keep_cpgs = NA_proportion_columns <= threshold_NA_cpg

  # Plot QC plots
  if (plot) {
    op = par(mfrow=c(2,2))

    hist(dnam$samples$missing, xlab="NA per sample", main=paste0("Discard ", sum(too_many_NAs_per_sample), "/", nrow(dnam$samples), " samples"))
    abline(v=threshold_NA_sample)

    hist(NA_proportion_columns, xlab="NA per CpG", main=paste0("Discard ", sum(!keep_cpgs, na.rm=T), "/", ncol(dnam$cpgs), " CpG sites"))
    abline(v=threshold_NA_cpg)

    boxplot(c(dnam$samples$bc1.grn, dnam$samples$bc1.red, dnam$samples$bc2) ~ rep(c("BC1 green", "BC1 red", "BC2"), each=nrow(dnam$samples)), ylim=c(0,50), main=paste0(sum(bc1_grn_failed, bc1_red_failed, bc2_failed), "/", nrow(dnam$samples), " samples failed bisulfite conversion"))
    abline(h=bc_threshold)

    plot(density(colMeans(eira$cpgs[keep_samples, keep_cpgs], na.rm=T), from=0, to=1), las=1, lwd=2, col="navy", main="Methylation β-value-distribution", xlab="")

    par(op)
  }

  # Actually remove samples and cpgs
  dnam$samples = dnam$samples[keep_samples, ]
  #dnam$excluded_cpgs = dnam$cpgs[rownames(dnam$samples), !keep_cpgs]
  dnam$cpgs = dnam$cpgs[rownames(dnam$samples), keep_cpgs]

  # Also remove corresponding samples from other list-entries in dnam (ratios, snps, possibly controls and intensities)
  dnam$snps = dnam$snps[rownames(dnam$samples), ]
  if (!is.null(dnam$intensities_A)) dnam$intensities_A = dnam$intensities_A[rownames(dnam$samples), ]
  if (!is.null(dnam$intensities_B)) dnam$intensities_B = dnam$intensities_B[rownames(dnam$samples), ]
  if (!is.null(dnam$controls_red)) dnam$controls_red = dnam$controls_red[rownames(dnam$samples), ]
  if (!is.null(dnam$controls_grn)) dnam$controls_grn = dnam$controls_grn[rownames(dnam$samples), ]

  dnam
}


infer_sex = function(dnam, kmeans=TRUE, centroids = rbind(m=c(0.3, 0.1), f=c(0.5, 0.7)), plot=TRUE) {
  # Determine sex by k-means with predefined cluster initialization
  if (kmeans) {
    # kmeans will throw an error when one cluster is empty, i.e. all samples are from the same sex
    sex = tryCatch({
      km = kmeans(cbind(dnam$samples$median.chrX, dnam$samples$missing.chrY), centroids)
      factor(km$cluster, levels=c(1,2), labels=c("m", "f"))
    }, error = function(err) {
      km = kmeans(cbind(dnam$samples$median.chrX, dnam$samples$missing.chrY), 1)
      if (dist(rbind(km$centers[1,], centroids["m",])) < dist(rbind(km$centers, centroids["f",]))) factor(rep("m", nrow(dnam$samples)), levels=c("m", "f"))
      else factor(rep("f", nrow(dnam$samples)), levels=c("m", "f"))
    })
    names(sex) = rownames(dnam$sample)

  # Alternatively determine sex by hard thresholds (advantage: classifies uncertain samples as NA)
  } else {
    sex = factor(rep(NA, times=nrow(dnam$samples)), levels=c("m", "f"))
    names(sex) = rownames(dnam$sample)

    threshold_chrX = mean(centroids[,1])
    threshold_chrY = mean(centroids[,2])

    # Female
    sex[rownames(dnam$samples)[dnam$samples$median.chrX > threshold_chrX &
                               dnam$samples$missing.chrY > threshold_chrY]] = "f"

    # Male
    sex[rownames(dnam$samples)[dnam$samples$median.chrX <= threshold_chrX &
                               dnam$samples$missing.chrY <= threshold_chrY]] = "m"
  }

  if (plot) {
    plot(dnam$samples$median.chrX, dnam$samples$missing.chrY, col=addNA(sex), xlim=c(0,1), ylim=c(0,1),
      xlab = "Median methylation level Chr. X",
      ylab = "Proportion missing Chr. Y")
    if (!kmeans) {
      rect(threshold_chrX, threshold_chrY, 1, 1, border = "red")
      rect(0.0, 0.0, threshold_chrX, threshold_chrY, border = "black")
    }
    legend("bottomright", c("male", "female"), col=c("black", "red"), pch=1)
  }

  sex
}


call_snps = function(dnam, non_carrier_threshold = 0.2, homozygous_threshold = 0.8, plot = TRUE) {
  calls = array(NA, dim=c(nrow(dnam$samples), ncol(dnam$snps)),
                dimnames=list(rownames(dnam$samples), colnames(dnam$snps)))

  calls[dnam$snps <= non_carrier_threshold] = 0 # non-carrier
  calls[dnam$snps > non_carrier_threshold & dnam$snps < homozygous_threshold] = 1 # heterozygous
  calls[dnam$snps >= homozygous_threshold] = 2 # homozygous

  if (plot) {
    title = paste(dim(dnam$snps)[1], "samples, ", dim(dnam$snps)[2], "SNPs")
    boxplot(as.numeric(as.matrix(dnam$snps)) ~ as.numeric(as.matrix(calls)), xlab="Carrier status", ylab="Theta intensities", main = title)
    abline(h=non_carrier_threshold)
    abline(h=homozygous_threshold)
  }

  calls
}

snp_distance = function(dnam, plot_histogram=TRUE, plot_heatmap=FALSE) {
  dist = dist(call_snps(dnam), method="manhattan")

  if (plot_histogram) {
    hist(dist, main="Genotype distance according to 65 SNPs", xlab="Absolute number of different alleles", freq=F)
  }

  if (plot_heatmap) {
    heatmap(as.matrix(dist), main="Genotype distance according to 65 SNPs")
  }

  dist
}

identify_replicates = function(dnam, snp_distance_margin = 0) {
  inferred.sex = infer_sex(dnam)
  snp.distance = snp_distance(dnam)

  # Identify possible replicate pairs: those with (near) zero SNP distance
  replicate_ids_in_lower_tri = which(snp.distance < snp_distance_margin)

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

  # Normalise so that rowSums are 1 and reflect the proportional nature of this data
  wbc.predictions * (1/rowSums(wbc.predictions))
}
