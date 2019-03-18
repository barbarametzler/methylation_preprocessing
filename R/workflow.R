# Load package
devtools::load_all()


### Preprocessing

# Preprocess all idat files in folder "idat"
eira = preprocess("hm450", "../idat", return.intensities=TRUE, verbose=TRUE)


### Quality control

# Remove samples with too many missing values (>10%) or failed bisulphate conversion
# Remove CpG sites with too many missing values (>20%)
eira_red = remove_unreliable_samples_probes(eira, threshold_NA_sample = 0.1, threshold_NA_cpg = 0.2, bc_outlier_SD = 3, plot=T)

# Infer sex
sex = infer_sex(eira, kmeans=T, plot=T)

# Call SNPs
snps = call_snps(eira, plot=T)
snp_distance(eira, plot_histogram=T, plot_heatmap=F)

# Identify replicates from the same individuals according to inferred sex and SNPs
identify_replicates(eira)

# Estimate leukocyte-subtype proportions
estimate_leukocytes(eira)
