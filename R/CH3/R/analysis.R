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
linear_mixed_models = function(dnam, outcome, chip, position, confounders=NA, lrt=TRUE, save.residuals=FALSE, save.ranks=TRUE, plot=FALSE) {
  # If dnam is a list as returned by this package's preprocess() function
  if (typeof(dnam) == "list") {
    chip = dnam$samples$chip
    position = dnam$samples$position
    dnam = dnam$cpgs
  }
  if (missing(chip) | missing(position)) {
    stop("Arguments chip and position have to be provided when dnam is a matrix and not the list output of preprocess()")
  }

  stopifnot(rownames(dnam) %in% names(outcome))

  data = data.frame(chip = chip, position = position, outcome = outcome[rownames(dnam)])
  formula_string = "dnam ~ outcome + (1|chip) + (1|position)"

  if (!is.na(confounders)) {
    stopifnot(rownames(dnam) %in% rownames(confounders))
    data = cbind(data, confounders[rownames(dnam),])
    formula_string = paste0(formula_string, "+", paste0(colnames(confounders), collapse="+"))
  }

  models = omics::mlmer(as.formula(formula_string), data, vars="outcome", lrt, save.residuals, save.ranks)

  if (plot) {
    omics::ranks.heatmap(models$ranef.ranks$chip)
    omics::ranks.heatmap(models$ranef.ranks$position)
  }

  models
}
