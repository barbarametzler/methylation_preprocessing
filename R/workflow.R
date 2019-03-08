# Load package
devtools::load_all()


### Preprocessing steps (according to Comp Epi Handbook)

# Preprocess all idat files in folder "idat"
eira = preprocess("hm450", "../idat", return.intensities=TRUE)


### Quality control steps (according to Comp Epi Handbook)

# Remove samples with too many missing values (>10%) and outliers
eira = remove_unreliable_samples(eira, thresholdNA = 0.1)

# Infer sex
infer_sex(eira, plot=T)

# Call SNPs
snps = call_snps(eira, plot=T)
snp_distance(eira, plot=T)

# Identify replicates from the same individual according to inferred sex and SNPs
identify_replicates(eira)

# Remove CpG-sites with too many missing values (>10%)
NA.proportion.columns = colMeans(eira$cpgs, na.rm=T)
table(NA.proportion.columns<0.1, useNA="always")
excluded_cpgs = eira$cpgs[, NA.proportion.columns >= 0.1]
eira$cpgs = eira$cpgs[, NA.proportion.columns < 0.1]
