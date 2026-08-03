# Create a memory effecient data structure for eQTL
# caculations.
slice_data = function(loci, expression_data) {
	load_it(c("MatrixEQTL", "lumi"))
	sliced_probes = SlicedData$new(as.matrix(loci))
	sliced_genes = SlicedData$new(as.matrix(expression_data))

	return(list(sliced_probes=sliced_probes, sliced_genes=sliced_genes))
} # end slice_data

# Perform the eQTL analysis using MatrixEQTL.
get_eqtls = function(loci, 
	expression_data, 
	probe_locs,
	gene_locs,
	sliced_dummy = SlicedData$new(as.matrix(character())),
	cis_cutoff = 1e-3,
	trans_cutoff = 1e-20,
	num_cis=1e+10,
	num_trans=50,
	use_cpg_islands=FALSE,
	orderby="beta") {

	if(use_cpg_islands) {
		by_y = "island"
	} else {
		by_y = "probe"
	}

	result = slice_data(loci, expression_data)
	sliced_probes = result$sliced_probes
	sliced_genes = result$sliced_genes

	result = eqtl(sliced_probes,
		sliced_genes,
		sliced_dummy,
		probe_locs,
		gene_locs,
		cis_cutoff,
		trans_cutoff
		)

	# result$cis$eqtls = result$cis$eqtls[result$cis$eqtls$FDR <= .05, ]
	# result$trans$eqtls = result$trans$eqtls[result$trans$eqtls$FDR <= .05, ]

	if(orderby=="beta"){
		result$cis$eqtls = result$cis$eqtls[with(result$cis$eqtls, order(-abs(beta))),]
		result$trans$eqtls = result$trans$eqtls[with(result$trans$eqtls, order(-abs(beta))),]
	}else{
		result$cis$eqtls = result$cis$eqtls[with(result$cis$eqtls, order(pvalue)),]
		result$trans$eqtls = result$trans$eqtls[with(result$trans$eqtls, order(pvalue)),]
	}

	result$cis$eqtls = result$cis$eqtls[!duplicated(result$cis$eqtls$snps),]

	if(nrow(result$cis$eqtls) > num_cis) {result$cis$eqtls=result$cis$eqtls[1:num_cis,]}

	cis_link = base::merge(result$cis$eqtls[,c(1,2)], probe_locs, by.x="snps", by.y=by_y)		
	cis_link = base::merge(cis_link, gene_locs, by.x="gene", by.y="gene")
	#cis_link = cis_link[match(eqtls$result$cis$eqtls$snps, cis_link$snps),]
	cislinks = paste0(cis_link[,2], cis_link[,1])
	ciseqtlids = paste0(result$cis$eqtls[,1], result$cis$eqtls[,2])
	cis_link = cis_link[match(ciseqtlids, cislinks),]
	
	result$trans$eqtls = result$trans$eqtls[!duplicated(result$trans$eqtls$gene),]

	if(nrow(result$trans$eqtls) > num_trans) {result$trans$eqtls=result$trans$eqtls[1:num_trans,]}

	trans_link = base::merge(result$trans$eqtls[,c(1,2)], probe_locs, by.x="snps", by.y=by_y)	
	trans_link = base::merge(trans_link, gene_locs, by.x="gene", by.y="gene")
	translinks = paste0(trans_link[,2], trans_link[,1])
	transeqtlids = paste0(result$trans$eqtls[,1], result$trans$eqtls[,2])
	trans_link = trans_link[translinks %in% transeqtlids,]

  	return(list(result=result, cis_link=cis_link, trans_link=trans_link))
} # end get eQTLs

eqtl = function(
	loci_values,
    trait_values,
    covar_values=NULL,
    loci_indices,
    trait_indices,
    cis,
    trans,
    cis_dist = 1e4) {

    'Performs eQTL analysis for arbitrary loci and traits (does not have to be SNPs).

    loci_values > some sort of value associated with each loci. This can be minor
      allele frequency for SNPs or methylation percent. Supply a data.frame with
      an id column with loci ids and then sample columns with values.
    trait_values > expression data in data.frame. First column should be id of gene
      followed by sample columns with expression values.
    covar_values > optional covariate data arranged in the same fashion
    loci_indices > first column is snp, gives SNP or other id, second column is chr
      for chromosome, given as "chr1" for example, third column is pos for position.
    trait_indices > geneid, chr, s1 (start), s2 (end)
    cis > significance threshold for cis eQTLs
    trans > significance threshold for trans eQTLs
    cis_dist > max cis distance to consider - has computational complexity consequences
    returns < MatrixEQTL object'
    
    load_it("MatrixEQTL")
	
	## Location of the package with the data files.
	base.dir = find.package('MatrixEQTL');
	
	## Settings
	
	# Linear model to use, modelANOVA, modelLINEAR, or modelLINEAR_CROSS
	useModel = modelLINEAR; # modelANOVA, modelLINEAR, or modelLINEAR_CROSS
	
	# Output file name
	output_file_name_cis = tempfile();
	output_file_name_tra = tempfile();
	
	# Only associations significant at this level will be saved
	pvOutputThreshold_cis = cis;
	pvOutputThreshold_tra = trans;
	
	# Error covariance matrix
	# Set to numeric() for identity.
	errorCovariance = numeric();
	
	# Distance for local gene-SNP pairs
	cisDist = cis_dist;
	
	snps = loci_values
	gene = trait_values
	cvrt = covar_values
	snpspos = loci_indices
	genepos = trait_indices
	
	me = Matrix_eQTL_main(
			snps = snps,
			gene = gene,
			cvrt = cvrt,
			output_file_name     = output_file_name_tra,
			pvOutputThreshold     = pvOutputThreshold_tra,
			useModel = useModel,
			errorCovariance = errorCovariance,
			verbose = TRUE,
			output_file_name.cis = output_file_name_cis,
			pvOutputThreshold.cis = pvOutputThreshold_cis,
			snpspos = snpspos,
			genepos = genepos,
			cisDist = cisDist,
			pvalue.hist = "qqplot",
			min.pv.by.genesnp = FALSE,
			noFDRsaveMemory = FALSE);
	
	unlink(output_file_name_tra);
	unlink(output_file_name_cis);
	
	## Results:
	return(me)
	#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
	#cat('Detected local eQTLs:', '\n');
	#show(me$cis$eqtls)
	#cat('Detected distant eQTLs:', '\n');
	#show(me$trans$eqtls)
	#return(me)
} #end eqtl

top_eqtl_genes = function(link, distance=1500) {
  locations = unique(paste0(link$chromosome, ":", link$location, "-", link$location))
  gene_info = query_genes_by_coordinates(locations, distance)
  return(gene_info)
} # end top_genes