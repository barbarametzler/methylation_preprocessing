# Load package
devtools::load_all()


### Preprocessing

# Preprocess all idat files in folder "idat"
eira = preprocess("hm450", "../idat", return_intensities=T, verbose=TRUE)
#saveRDS(eira, "eira_list.rds")

### Quality control

# Remove samples with too many missing values (>10%) or failed bisulfite conversion
# Remove CpG sites with too many missing values (>10%)
eira_red = remove_unreliable_samples_probes(eira, threshold_NA_sample = 0.1, threshold_NA_cpg = 0.1, bc_threshold = 1, plot=T)

# Infer sex
inferred_sex = infer_sex(eira, kmeans=T, plot=T)

# Call SNPs
snps = call_snps(eira, plot=T)
system.time(snp_distance(eira, plot_histogram=T, plot_heatmap=F))

# Identify replicates from the same individuals according to inferred sex and SNPs
identify_replicates(eira)

# Estimate leukocyte-subtype proportions
estimate_leukocytes(eira)



### Analysis

# Read covariates and table containing ids & set rownames of covariates to match id-format of dnam-dataframes
covars = readRDS("../Covariates.rds")
sample_sheet = readRDS("../Sample_sheet.rds")
rownames(covars) = sample_sheet$sample.id

# Compare inferred sex with prior knowledge
supplied_sex = covars$gender
names(supplied_sex) = rownames(covars)
compare_sex(eira, supplied_sex)

# Linear mixed models
outcome = covars$disease_status == "rheumatoid arthritis"
names(outcome) = rownames(covars)
models = linear_mixed_models(eira, outcome=outcome, confounders=covars[,3:5], plot=F)
saveRDS(models, "models.rds")

volcano_plot(models$coefficients, T)
manhattan_plot(models$coefficients, "hm450", T)
