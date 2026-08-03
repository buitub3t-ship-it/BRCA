#impute.R
library(impute)

impute_meth = function(df, na_max_prop=.5) {
	load_it("lumi")

	# Filter probes more than half NAs
	df = df[rowMeans(is.na(df)) <= na_max_prop,]
	#df = beta2m(df)

	df = data.frame(impute.knn(as.matrix(df), k=10, rng.seed=3301)$data)

	return(df)
} # end impute_meth

