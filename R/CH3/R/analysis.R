#' Performs multiple linear mixed models
#'
#' One model is fitted per CpG site. Beta-values are treated as dependent variables and outcomes as independent variables. Chip id and position are treated as random effect intercepts in the mixed model to reduce technical variance. Additional confounders like age and sex can optionally be added and treated as fixed effects.
#'
#' @param dnam list-output of \link{preprocess} or generic methylation-beta-value-matrix
#' @param outcome named vector of outcome variable
#' @param chip only has to be provided when \code{dnam} is not a list-output of \link{preprocess}. Vector of chip id each sample was performed on.
#' @param position only has to be provided when \code{dnam} is not a list-output of \link{preprocess}. Vector of chip position each sample was performed on.
#' @param confounders optional dataframe with confounders, rownames matching those of dnam-list-entries. Even when only one confounder is present provide it as a 1-column-dataframe.
#' @param lrt whether p-values should be computed using likelihood-ratio tests. See ?omics::mlmer
#' @param save.residuals whether model residuals should be stored and returned. See ?omics::mlmer
#' @param save.ranks whether random effect ranks should be stored and returned. See ?omics::mlmer
#' @param plot \code{TRUE} plots two heatmap for the random effect ranks of chip id and position. See ?omics::ranks.heatmap
#'
#' @return A list with elements \code{coefficients} (containing \code{coef}, \code{coef.se} and \code{pval} for each CpG site) and optionally \code{residuals} and \code{ranef.ranks}. See ?omics::mlmer
#'
#' @seealso \link{preprocess}
#'
#' @export
linear_mixed_models = function(dnam, outcome=NA, confounders=NA, chip, position, lrt=TRUE, save.residuals=FALSE, save.ranks=TRUE, plot=FALSE) {
  # If dnam is a list as returned by this package's preprocess() function
  if (typeof(dnam) == "list") {
    chip = dnam$samples$chip
    position = dnam$samples$position
    dnam = dnam$cpgs
  }
  if (missing(chip) | missing(position)) {
    stop("Arguments chip and position have to be provided when dnam is a matrix and not the list output of preprocess()")
  }

  data = data.frame(position = position)
  formula_string = "dnam ~ (1|position)"

  # Data from more than one chip?
  if (length(levels(as.factor(chip))) > 1) {
    data = cbind(data, chip)
    formula_string = paste(formula_string, "+ (1|position)")
  }

  if (!is.na(confounders)) {
    stopifnot(rownames(dnam) %in% rownames(confounders))
    data = cbind(data, confounders[rownames(dnam),])
    formula_string = paste(formula_string, "+", paste0(colnames(confounders), collapse="+"))
  }

  if (!is.na(outcome)) {
    stopifnot(rownames(dnam) %in% names(outcome))
    data = cbind(data, outcome=outcome[rownames(dnam)])
    formula_string = paste(formula_string, "+ outcome")

    models = omics::mlmer(as.formula(formula_string), data, vars="outcome", lrt, save.residuals, save.ranks)
  } else {
    # When no outcome is specified use LMM for denoising and return residuals no matter what provided argument
    models = omics::mlmer(as.formula(formula_string), data, lrt=FALSE, save.residuals=TRUE, save.ranks=save.ranks)
  }

  if (plot) {
    omics::ranks.heatmap(models$ranef.ranks$chip)
    omics::ranks.heatmap(models$ranef.ranks$position)
  }

  models
}

volcano_plot = function(results, annotate = FALSE, threshold = 0.05) {
  fdr = p.adjust(results$pval, method = "BH")

  par(mar = c(4.5, 4.5, 1, 1))
  plot(results$coef, -log10(results$pval), pch = 19,
       las = 1, cex = 0.5, xlab = expression(beta),
       ylab = expression(-log[10](p[value])), col = ifelse(fdr < threshold, yes = "tomato", no = "darkgrey"))
  abline(v = 0, lty = 3)
  abline(h = -log10(threshold/nrow(results)), lty = 2, col = "darkred")
  legend("bottomleft", col = c("darkred", "tomato", "darkgrey"), lty = c(2, NA, NA), pch = c(NA, 19, 19), cex = 0.7,
         legend = c("Bonferroni threshold", "FDR significant hits", paste("Not significant at", threshold)))
  if (annotate) {
    text(results$coef, -log10(results$pval), pos = 3,
         offset = 0.2, cex = 0.5, labels = ifelse(fdr < threshold, yes = rownames(results), no = ""))
  }
}

manhattan_plot = function(results, platform = NA, annotate = FALSE, threshold = 0.05) {
  stopifnot(platform %in% c("hm450", "epic"))
  probes=readRDS(paste0('R/illumina_manifests/', sprintf("%s_probes.rds", platform)))

  colors = c("#A6CEE3", "#65A4CC", "#247BB6", "#5EA4A1", "#A5D68D", "#80C665", "#43A838", "#789D51",
             "#D89B86", "#F37372", "#E83537", "#E94531", "#F69359", "#FDB156", "#FE9221", "#F58725",
             "#DCA08B", "#BDA2CE", "#8F6AB1", "#764D99", "#BEAA99", "#FBF794", "#D6A85E", "#B15928")

  par(mar = c(3, 4.5, 1, 1))
  plot(probes[rownames(results), "pos"], -log10(results$pval),
       col = colors[as.integer(probes[rownames(results), "chr"])], pch = 19, cex = 0.5,
       xlab = "", ylab = expression(-log[10](p[value])),
       las = 1, xaxt = "n")
  #axis(side = 1, at = chr_boundaries[-length(chr_boundaries)] +
  #       (chr_boundaries[-1] - chr_boundaries[-length(chr_boundaries)])/2,
  #     labels = c(seq(1:22), "X", "Y")[1:24], tick = FALSE,
  #     cex.axis = 0.6)
  #axis(side = 1, at = chr_boundaries, labels = NA,
  #     col.ticks = "grey")
  abline(h = -log10(threshold/nrow(results)), lty = 2,
         col = "darkred")
  if (annotate) {
    text(probes[rownames(results), "pos"], -log10(results$pval), pos = 3, offset = 0.2, cex = 0.5,
         labels = ifelse(results$pval < threshold/nrow(results), yes = rownames(results), no = ""))
  }
}
