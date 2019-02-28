#e1071::cmeans
#mvoutlier::pcout

remove.unreliable.samples = function(dnam, thresholdNA = 0.1) {
  # Remove samples with too many missing values
  dnam$samples = dnam$samples[dnam$samples$missing < thresholdNA, ]

  # Remove potential multivariate outliers according to BC controls
  tech.vars = scale(dnam$samples[,c("bc1.grn", "bc1.red", "bc2")])
  outliers <- mvoutlier::pcout(tech.vars)
  dnam$samples = dnam$samples[outliers$wfinal01 == 1, ]

  # Also remove corresponding samples from other list-entries in dnam (ratios, snps, possibly controls and intensities)
  dnam$cpgs = dnam$cpgs[dnam$samples$sample.id, ]
  dnam$snps = dnam$snps[, dnam$samples$sample.id, ]
  if (!is.null(dnam$intensities)) dnam$intensities = dnam$intensities[, dnam$samples$sample.id, ]
  if (!is.null(dnam$controls)) dnam$controls = dnam$controls[, dnam$samples$sample.id, ]

  dnam
}


infer.sex = function(dnam) {
  tmp <- e1071::cmeans(dnam$samples[,c("median.chrX", "missing.chrY")], centers=2)
  membership <- apply(tmp$membership, 1, function(x) {
    if (any(x > 0.95)) {
      which.max(x)
    } else {
      NA
    }
  })
  idx.m <- which.min(tmp$centers[,"missing.chrY"])
  idx.f <- which.max(tmp$centers[,"missing.chrY"])

  # Return inferred sex
  factor(membership, c(idx.m, idx.f), c("M", "F"))
}


call.snps = function(dnam) {
  membership = lapply(apply(dnam$snps["theta", dimnames(dnam$snps)[[2]],], 2, function(x) {
    idx <- !is.na(x)
    tmp <- e1071::cmeans(x[idx], centers=c(0, 0.5, 1))
    results <- matrix(NA, length(x), 3, dimnames=list(names(x), 0:2))
    results[which(idx),] <- tmp$membership
    results[which(!idx),] <- 1/3
    list(results)
  }), `[[`, 1)

  calls = do.call(cbind, lapply(membership, function(x) {
    apply(x, 1, function(x) {
      if (any(x > 0.95)) {
        which.max(x) - 1
      } else {
        NA
      }
    })
  }))

  calls
}


identify.replicates = function(dnam) {
  samples = dnam$samples
  samples$inferred.sex = infer.sex(dnam)
  calls = call.snps(dnam)

  # Merge SNPs
  samples <- cbind(samples, calls)

  # Compute similarity matrix
  K <- matrix(0, nrow(samples), nrow(samples),
              dimnames=list(rownames(samples), rownames(samples)))

  tmp <- t(data.matrix(samples[,c("inferred.sex", colnames(calls))]))

  for (i in 1:nrow(K)) {
    K[i,] <- as.integer(apply((tmp[,i] - tmp) == 0, 2, all, na.rm=TRUE))
  }
  table(K)
  table(K[lower.tri(K)])

  # Identify duplicates from SNPs
  dups <- as.data.frame(which(K == 1, arr.ind=TRUE))
  dups$sample1 <- rownames(samples)[dups$row]
  dups$sample2 <- rownames(samples)[dups$col]
  dups$id1 <- samples[dups$row,"sample.id"]
  dups$id2 <- samples[dups$col,"sample.id"]
  dups$source <- "snps"
  dups$row = dups$col <- NULL

  dups <- subset(dups,
                 sample1 != sample2
  )

  dups
}
