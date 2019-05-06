#quadprog::solve.QP


#' Removes unreliable samples and CpG sites
#'
#' Samples are deemed unreliable if they have too many NA values in CpG matrix or have failed bisulfite conversion according to 3 summary statistics in samples dataframe (bc1.grn, bc1.red, bc2). CpG sites are deemed unreliable if they have too many NA values.
#'
#' @param dnam list-output of \link{preprocess}
#' @param threshold_NA_sample samples (= rows) with more NAs than this threshold-fraction in their CpG matrix will be removed (default 0.1)
#' @param threshold_NA_cpg CpG sites (= columns) with more NAs than this threshold-fraction will be removed (default 0.2)
#' @param bc_threshold samples that failed bisulfite conversion (bc) according to this lower-threshold for bc-summary-statistics in samples dataframe will be removed (default 1)
#' @param plot \code{TRUE} for quality control plots showing NA fractions per row/column, bc-failure and average methylation beta-value-distribution among all samples
#'
#' @return Returns list in same format as \link{preprocess}
#'
#' @references Illumina BeadArray Controls Reporter Software Guide (Document #1000000004009 v00, October 2015): https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf
#'
#' @seealso \link{preprocess}
#'
#' @export
remove_unreliable_samples_probes = function(dnam, threshold_NA_sample = 0.1, threshold_NA_cpg = 0.2, bc_threshold = 1, plot=T) {
  # Define samples with failed bisulfite conversion (Illumina recommends threshold of 1: Illumina BeadArray Controls Reporter Software Guide (Document #1000000004009 v00, October 2015): https://support.illumina.com/content/dam/illumina-support/documents/documentation/chemistry_documentation/infinium_assays/infinium_hd_methylation/beadarray-controls-reporter-user-guide-1000000004009-00.pdf)
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

    plot(density(colMeans(dnam$cpgs[keep_samples, keep_cpgs], na.rm=T), from=0, to=1), las=1, lwd=2, col="navy", main="Avg. methylation β-value-distribution", xlab="")

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

kmeans_default = TRUE
centroids_default = rbind(m=c(0.3, 0.1), f=c(0.5, 0.7))

#' Infer sex
#'
#' Sex inference is performed based on two summary statistics: and median.chrX (higher for females because of X-chromosome-inactivation) and missing.chrY (higher for females but not 100% because of cross-hybridization of X-chromosomal-genes).
#'
#' @param dnam list-output of \link{preprocess}
#' @param kmeans \code{TRUE} performs k-means-clustering with predefined initialization. \code{FALSE} uses hard boundaries (\code{colMeans(centroids)}) which also yields NA values when uncertain.
#' @param centroids default initialization for k-means-clustering and boundaries
#' @param plot \code{TRUE} shows scatter plot with sex-inference
#'
#' @return Returns named factor-vector with levels "m" / "f" / NA (for \code{kmeans=FALSE} only)
#'
#' @seealso \link{preprocess}, \link{compare_sex}
#'
#' @export
infer_sex = function(dnam, kmeans=kmeans_default, centroids=centroids_default, plot=TRUE) {
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

#' Compares inferred sex with specified sex
#'
#' @param dnam list-output of \link{preprocess}
#' @param specified_sex named vector with prespecified sex
#' @param kmeans passed on to \link{infer_sex}
#' @param centroids passed on to \link{infer_sex}
#' @param plot \code{TRUE} shows scatter plot with specified sex
#'
#' @return Returns confusion matrix showing accordance of inferred and specified sex
#'
#' @seealso \link{preprocess}, \link{infer_sex}
#'
#' @export
compare_sex = function(dnam, specified_sex, kmeans=kmeans_default, centroids=centroids_default, plot=T) {
  stopifnot(rownames(dnam$samples) %in% names(specified_sex))

  specified_sex = specified_sex[rownames(dnam$samples)]
  inferred_sex = infer_sex(dnam, kmeans, centroids, plot=F)

  if (plot) {
    plot(dnam$samples$median.chrX, dnam$samples$missing.chrY, col=addNA(specified_sex), xlim=c(0,1), ylim=c(0,1),
         xlab = "Median methylation level Chr. X",
         ylab = "Proportion missing Chr. Y")
    if (!kmeans) {
      rect(threshold_chrX, threshold_chrY, 1, 1, border = "red")
      rect(0.0, 0.0, threshold_chrX, threshold_chrY, border = "black")
    }
    legend("bottomright", levels(specified_sex), col=c("black", "red"), pch=1)
  }

  table(inferred_sex, specified_sex)
}

#' Call Single Nucleotide Polymorphisms (SNPs)
#'
#' Even though the Illumina Infinium Bead Chips are made for CpG sites, they also contain a few (~60) highly polymorphic SNPs to help identify replicates. This function "calls" these SNPs, meaning it defines them as either non-carrier (0), heterozygous (1) or homozygous (2) based on the theta-value in the polar coordinate system.
#'
#' @param dnam list-output of \link{preprocess}
#' @param non_carrier_threshold SNPs with theta below this will be classified as non-carrier (default 0.2)
#' @param homozygous_threshold SNPs with theta above this will be classified as homozygous (default 0.8)
#' @param plot \code{TRUE} shows box plot-overview of theta-values among the three classes for all SNPs
#'
#' @return Returns matrix with samples in rows and SNPs in columns; 0 for non-carrier, 1 for heterozygous, 2 for homozygous
#'
#' @seealso \link{preprocess}, \link{snp_distance}
#'
#' @export
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

#' Computes genotype distance
#'
#' L1-distance-matrix of samples according to the ~60 SNPs present on Illumina Infinium Bead Chips corresponds to absolute number of different alleles.
#'
#' @param dnam list-output of \link{preprocess}
#' @param plot_histogram \code{TRUE} shows distribution of genotype-distance
#' @param plot_heatmap \code{TRUE} shows heatmap with all samples
#'
#' @return Returns \code{dist}-object which is lower triangle of distance matrix. To convert to matrix call \code{as.matrix()}
#'
#' @seealso \link{preprocess}, \link{call_snps}
#'
#' @export
snp_distance = function(dnam, plot_histogram=TRUE, plot_heatmap=FALSE) {
  dist = dist(call_snps(dnam, plot=F), method="manhattan")

  if (plot_histogram) {
    hist(dist, main="Genotype distance according to all SNPs", xlab="Absolute number of different alleles", freq=F)
  }

  if (plot_heatmap) {
    heatmap(as.matrix(dist), main="Genotype distance according to all SNPs")
  }

  dist
}

#' Identifies replicates
#'
#' Two samples are considered replicates when their inferred sex is the same and their SNP-L1-distance doesn't exceed \code{snp_distance_margin}.
#'
#' @param dnam list-output of \link{preprocess}
#' @param snp_distance_margin specifies the maximum number of different alleles for two samples to still be considered replicate (default 0)
#' @param kmeans passed on to \link{infer_sex}
#' @param centroids passed on to \link{infer_sex}
#'
#' @return Returns dataframe with replicates
#'
#' @seealso \link{preprocess}, \link{call_snps}, \link{snp_distance}, \link{infer_sex}
#'
#' @export
identify_replicates = function(dnam, snp_distance_margin = 0, kmeans=kmeans_default, centroids=centroids_default) {
  inferred.sex = infer_sex(dnam, kmeans, centroids, plot=F)
  snp.distance = snp_distance(dnam, plot_histogram=F)

  # Identify possible replicate pairs: those with (near) zero SNP distance
  replicate_ids_in_lower_tri = which(snp.distance <= snp_distance_margin)

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

#' Estimate leukocyte proportions
#'
#' Estimates white-blood-cell proportions from blood samples using the Houseman method and coefficients derived from the Reinius dataset. The following subtypes are estimated: Neutrophils, eosinophils, mast cells, CD4 T-cells, CD8 T-cells, B-cells, NK-cells
#'
#' @param dnam list-output of \link{preprocess}
#'
#' @return Returns a dataframe with samples in rows and leukocyte-subtype-proportions in columns.
#'
#' @references Houseman, Eugene Andres, et al. "DNA Methylation Arrays as Surrogate Measures of Cell Mixture Distribution." BMC Bioinformatics, vol. 13, no. 1, 2012, doi:10.1186/1471-2105-13-86.
#'
#' Reinius, Lovisa E., et al. "Differential DNA Methylation in Purified Human Blood Cells: Implications for Cell Lineage and Studies on Disease Susceptibility." PLoS ONE, vol. 7, no. 7, 2012, doi:10.1371/journal.pone.0041361.
#'
#' @seealso \link{preprocess}
#'
#' @export
estimate_leukocytes = function(dnam) {
  coefs = readRDS('R/illumina_manifests/leukocyte_coefs.rds')

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
