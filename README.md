## M2EFM

The M2EFM package allows you to build prognostic models derived from
gene expression features that are associated with disrupted methylation
patterns and clinical features. The final model is a Cox proportional
hazards model with a molecular risk score labelled "pred" in the model
and clinical covariates, allowing for easy interpretation.

To build an M2EFM model, you minimally need a file with DNA methylation
Beta-values from an Illumina 450k assay, a file with gene expression
values, and a file listing a set of methylation probes to use. The probes
should be sorted into some order that is meaningful to your experiment
(such as by p-value).

```{r, echo=TRUE, results="hide"}
	library(M2EFM)
	# The following lines take a path to a file (here provided to files in the 
	# package) and return MethylationData and ExpressionData objects, which
	# are used by some of the functions in this package.
	meth = load_meth(system.file("extdata", "TCGA_BRCA_METH.txt", package="M2EFM"))
	exp = load_exp(system.file("extdata", "TCGA_BRCA_EXP.txt", package="M2EFM"))

	# Get m2eQTLs.
	m2e = get_m2eqtls(probe_list=system.file("extdata", "PROBELIST.txt", package="M2EFM"), 
		meth_data=meth, 
		exp_data=exp)

	# Load the main clinical data.
	clin = read.table(system.file("extdata", "TCGA_BRCA_CLIN.txt", package="M2EFM"), 
		sep="\t", header=T, row.names=1)

	# Set covariates to be used in models.
	covariates = c("pathologic_stage", "age.Dx")

	# Assuming this is TCGA data, filter out all but the non-metastatic tumor samples.
	exp$values = exp$values[,grep("-01", colnames(exp$values))]
	clin = clin[grep("-01", rownames(clin)),]

	# Get genes and probes from m2eQTL models.
	genes = get_genes(m2e, "trans")
	exp$values = exp$values[genes,]

	probes = get_probes(m2e)
	meth$m_values = meth$m_values[probes,]

	# Create combined dataset of genes and probes.
	comb = merge(t(exp$values), t(meth$m_values), by=c("row.names"))
	rownames(comb) = comb[,1]
	comb = comb[,-1]

	# Subset clinical data to same samples as combined data.
	clin = clin[rownames(comb),]

	# Create a survival object.
	surv = Surv(clin$OVERALL.SURVIVAL, clin$overall.survival.indicator)

	# Only keep samples with outcomes.
	has_surv = !is.na(surv)
	surv = surv[has_surv]
	clin = clin[has_surv,]
	comb = comb[has_surv,]

	# Build M2EFM model. This function takes a data.frame or matrix, because
	# the data may not be just expression or methylation (it may be both).
	m2efm = build_m2efm(comb, surv, clin, covariates, standardize=TRUE)
```

```{r, echo=TRUE}
	summary(m2efm)
```

## Validation

The M2EFM package also allows you to evaluate your model using internal and
external validation data, using the evaluate() function. For internal validation,
the data are split into random training and test datasets, with the number of
random splits determined by monte_carlo_size. The proportion of the split is
set by training_prop. If the validation data do not have both methylation and
gene expression, set integrate_data to FALSE to prevent a second data integration.

```{r, echo=TRUE, results="hide", message=FALSE, warning=FALSE}
	library(M2EFM)
	meth = load_meth(system.file("extdata", "TCGA_BRCA_METH.txt", package="M2EFM"))
	exp = load_exp(system.file("extdata", "TCGA_BRCA_EXP.txt", package="M2EFM"))
	m2e = get_m2eqtls(probe_list=system.file("extdata", "PROBELIST.txt", package="M2EFM"), 
		meth_data=meth, 
		exp_data=exp)

	# This time we are going to use some external validation data.
	kao = load_exp(system.file("extdata", "KAO_BRCA_EXP.txt", package="M2EFM"))

	clin = read.table(system.file("extdata", "TCGA_BRCA_CLIN.txt", package="M2EFM"), 
		sep="\t", header=T, row.names=1)

	kao_clin = read.table(system.file("extdata", "KAO_BRCA_CLIN.txt", package="M2EFM"), 
		sep="\t", header=T, row.names=1)

	covariates = c("pathologic_stage", "age.Dx")

	exp$values = exp$values[,grep("-01", colnames(exp$values))]
	clin = clin[grep("-01", rownames(clin)),]

	# Build and evaluate the M2EFM model.
	eval_results = evaluate(
		meth, exp, clin = clin, eqtls = m2e,  
		gene_type="trans", gene_sig=GENE_SIG, 
		covariates=covariates, 
		monte_carlo_size=100, 
		training_prop = .7,
		event_col = "overall.survival.indicator",
		time_to_event_col = "OVERALL.SURVIVAL",
		integrate_data=FALSE,
		val=list(kao$values),
		val_clin=list(kao_clin))

	plot(eval_results)
	summary(eval_results)
```
# M2EFM
