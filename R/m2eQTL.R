
# source("dna_methylation.R")
# source("diff_meth.R")
# source("MethylationData.R")
# source("ProfileData.R")
# source("genome.R")
# source("ExpressionData.R")
# source("eqtl.R")
# source("impute.R")

#' Finds m2eQTLs from gene expression and DNA methylation data.
#'
#' @param probe_list_loc path to a list of probe ids, used to filter methylation data
#' @param meth_data path to a matrix of methylation data or a MethylationData object
#' @param exp_data path to expression data or an ExpressionData object
#' @param gene_list a file with a list of genes to filter to (optional)
#' @param num_cis number of cis m2eQTLs to return
#' @param num_trans number of trans m2eQTLs to return
#' @param num_probes number of top probes to use in analysis
#' @export
get_m2eqtls = function(
	probe_list_loc = NULL,
	meth_data = NULL,
	exp_data = NULL,
	gene_list = NULL,
	num_cis = 1e+10,
	num_trans = 100,
	num_probes = 300,
	treatment_col = "-01",
	cc=1e-3,
	tc=1e-10
) {
	# Check for required parameters.
	if(is.null(probe_list_loc) || is.null(meth_data) || is.null(exp_data)) {
		stop("probe_list_loc, meth_data, and exp_data must be provided")
	}
	
	# Get the list of probes differentially methylated.
	probe_list = read.table(probe_list_loc, header=F, row.names=NULL)
	probe_list = probe_list[,1]

	if("data.frame" %in% class(meth_data)) {
		meth_data = MethylationData$new(meth_data)
	} else if(!"MethylationData" %in% class(meth_data)) {
		meth_data = load_meth(meth_data)
	} 

	# Remove from the probe list any probes that are not in the data (there
	# shouldn't be many at this point, just those with at least 50% missing
	# data).
	probe_list = probe_list[as.character(probe_list) %in% rownames(
		meth_data$beta_values)]

	# Subset the data to those probes left in the probe list.
	meth_data$beta_values = meth_data$beta_values[as.character(probe_list),]

	# Get a table of locations for the probes.
	probe_locs = get_probe_locs(as.character(probe_list))
	probe_locs$probe = probe_list

	if("data.frame" %in% class(exp_data)) {
		exp_data = ExpressionData$new(exp_data)
	} else if(!"ExpressionData" %in% class(exp_data)) {
		exp_data = load_exp(exp_data)
	}

	# Keep only tumor samples.
	exp_data$values = exp_data$values[,grep(treatment_col, colnames(exp_data$values))]

	# Filter out unexpressed genes.
	exp_cancer = mad_filter(exp_data, thres=.05)$values

	if(substr(rownames(exp_cancer)[1], 1, 4) == "ENST") {
		genes = readRDS(system.file("extdata", "EnsemblTranscripts.RDS", package="M2EFM"))
		transcripts = substr(rownames(exp_cancer), 1, 15)
		genes = genes[!duplicated(genes[,1]),]
		rownames(genes) = genes[,1]
		exp_cancer = exp_cancer[transcripts %in% rownames(genes),]
		transcripts = substr(rownames(exp_cancer), 1, 15)
		genes = genes[,c(1,3,4,5)]
		colnames(genes) = c("gene", "chrom", "start", "end")
		genes$chrom = paste0("chr", genes$chrom)
		head(genes[transcripts,])
		gene_locs = genes[transcripts,]
	
	} else {
		entrez_ids = get_entrez_ids(rownames(exp_cancer))
		
		# Subset expression data to those with entrez ids.
		exp_cancer = exp_cancer[as.character(entrez_ids[,1]),]
		
		# Get table of gene locations.
		gene_locs = get_gene_locs(as.character(entrez_ids[,2]))
		gene_locs = unique(gene_locs)
		
		complete = entrez_ids[entrez_ids[,1] %in% gene_locs[,1],]
		
		gene_locs = gene_locs[gene_locs$gene %in% rownames(exp_cancer),]
		
		exp_cancer = exp_cancer[gene_locs$gene,]
	
	}


	# Filter methylation data down to those with differentially
	# methylated probes from the discovery set, which were passed 
	# to this routine.
	diff_probes = meth_data$beta_values[as.character(probe_list),grep(treatment_col, colnames(meth_data$beta_values))]

	colnames(diff_probes) = substr(colnames(diff_probes),1,15)

	#diff_probes = meth_data$beta_values[as.character(probe_list),]
	diff_probes = diff_probes[,colnames(diff_probes) %in% colnames(exp_cancer)]

	# Filter expression data to those samples with methylation data.
	exp_cancer = exp_cancer[,as.character(colnames(diff_probes))]

	probe_data = diff_probes[1:num_probes,]

	# normals = grep("-11", colnames(probe_data))
	# tns = probe_data[,substr(colnames(probe_data),1,12) %in% substr(colnames(probe_data)[normals],1,12)]

	# Find m2eQTLs.
	eqtls = get_eqtls(as.matrix(beta2m(probe_data)), 
		as.matrix(cbind(exp_cancer)), 
		probe_locs, 
		gene_locs, 
		num_cis=num_cis, 
		num_trans=num_trans, 
		use_cpg_islands=FALSE, 
		orderby="beta",
		cis_cutoff=cc,
		trans_cutoff=tc)

	#eqtls = get_eqtls(as.matrix(beta2m(tns)), as.matrix(cbind(exp_cancer[,as.character(colnames(tns))])), probe_locs, gene_locs, num_cis=700, num_trans=700, use_cpg_islands=FALSE, orderby="beta")

	#return(list(cis=eqtls$result$cis$eqtls, trans=eqtls$result$trans$eqtls))
	return(eqtls)
} # end perform_workflow

#' Takes a vector of methylation probe ids and builds a table of genomic 
#' probe loci. 
#'
#' @param probe_list vector of probe ids
#' @export
get_probe_locs = function(probe_list) {
	# # Query the 450k annotation for probe locations.
	# probe_locs = query_450k_annotation(
	# 	columns=c("IlmnID", "CHR", "MAPINFO"), 
	# 	probe_ids=probe_list)
	# colnames(probe_locs) = c("probe", "chromosome", "location")

	#load_it("IlluminaHumanMethylation450kanno.ilmn12.hg19")

	locs = getLocations(IlluminaHumanMethylation450kanno.ilmn12.hg19)

	probe_locs = as.data.frame(locs[as.character(probe_list),])
	
	probe_locs = data.frame(probe=rownames(probe_locs),
		chromosome=probe_locs$seqnames, location=probe_locs$start)

	# # Fix the chromosome labels.
	# probe_locs$chromosome = paste0("chr", probe_locs$chromosome)

	return(probe_locs)
} # end get_probe_locs()

#' Takes a vector of gene symbols and converts them to entrez ids.
#'
#' @param gene_symboles a vector of gene symbols
#' @return a data.frame of gene symbols and corresponding entrez ids
#' @export
get_entrez_ids = function(gene_symbols) {
	load_it("org.Hs.eg.db")

	# Find gene ids for genes by their symbols.
	x = org.Hs.eg.db::org.Hs.egALIAS2EG
	mapped_genes = mappedkeys(x)
	xx = AnnotationDbi::as.list(x[mapped_genes])

	found_ids = gene_symbols %in% names(xx)
	entrez_ids = sapply(xx[gene_symbols[found_ids]], function(x) x[1])

	return(data.frame(gene=gene_symbols[found_ids], entrez=entrez_ids))
} # end get_entrez_ids()

#' Takes a data.frame with entrez ids and gene symbols and builds a table of 
#' gene locations.
#' @param entrez_ids data.frame with first column gene symbols and second column the corresponding entrez ids
#' @return a data.frame with gene locations 
#' @export
get_gene_locs = function(entrez_ids) {
	require(org.Hs.eg.db)
	require(AnnotationDbi)

	# Get symbol map
	sym_map = org.Hs.egSYMBOL
    # Get the gene symbol that are mapped to an entrez gene identifiers
    mapped_genes = mappedkeys(sym_map)
    # Convert to a list
    gene_list = AnnotationDbi::as.list(sym_map[mapped_genes])

	entrez_ids = entrez_ids[entrez_ids %in% names(gene_list)]
	if(length(entrez_ids) == 0) return(data.frame()) 
	gene_list = gene_list[entrez_ids]

	# Get TSS map
    loc_map = org.Hs.egCHRLOC
    # Get the entrez gene identifiers that are mapped to chromosome locations
    mapped_genes = mappedkeys(loc_map)
    # Convert to a list
    tss_list = AnnotationDbi::as.list(loc_map[mapped_genes])
	tss_strand = sapply(tss_list[entrez_ids], function(x) x>0)
	has_na = sapply(tss_strand, function(x) x[1])
	has_na = is.na(has_na)
	tss_strand = tss_strand[!has_na]

	tss_list = sapply(tss_list, abs)

	tss_list = tss_list[entrez_ids]

	gene_list = gene_list[!has_na]
	tss_list = tss_list[!has_na]

	# Get gene end map
	end_map = org.Hs.egCHRLOCEND
	# Get the entrez gene identifiers that are mapped to chromosome locations
	mapped_genes = mappedkeys(end_map)
	# Convert to a list
    end_list = AnnotationDbi::as.list(end_map[mapped_genes])
	end_list = sapply(end_list, abs)

	end_list = end_list[entrez_ids]

	end_list = end_list[!has_na]

	start_list = sapply(1:length(tss_list), function(i) {
		if(tss_strand[[i]][1]){
			min(tss_list[[i]])
		} else {
			max(tss_list[[i]])
		}
	})

	end_list = sapply(1:length(end_list), function(i) {
		if(tss_strand[[i]][1]){
			max(end_list[[i]])
		} else {
			min(end_list[[i]])
		}
	})

	chrom_list = sapply(tss_list, function(x) names(x)[1])
	entrez_ids = names(chrom_list)
	chromosomes = paste0("chr", unlist(chrom_list))

	symbols = unlist(gene_list[entrez_ids])

	tss_strands = sapply(tss_strand, function(x) x[1])

	gene_locs = data.frame(gene=symbols, chrom=chromosomes, start=start_list, end=end_list, pos_strand=unlist(tss_strands))

	for(i in 1:nrow(gene_locs)) {
		if(!gene_locs[i,]$pos_strand) {
			tmp = gene_locs[i,]$start
			gene_locs[i,]$start = gene_locs[i,]$end
			gene_locs[i,]$end = tmp
		}
		if(is.na(gene_locs[i,]$gene)) {
			gene_locs[i,]$gene = entrez_ids[i]
		}
	}

	# gene_locs = query_gene_coordinates(entrez_ids[,2], min=T, max=T)

	# # Add gene names to gene locations.
	# gene_locs = cbind(gene_locs, gene=entrez_ids[,1])

	# # This is a hack to fix coordinates from UCSC.
	# gene_locs[gene_locs$gene=="TNXB",5] = 32064150
	# gene_locs[gene_locs$gene=="KIFC1",5] = 33359763
	# gene_locs[gene_locs$gene=="CRHR1",5] = 43907551

	# # Cut down to only the four columns needed for eqtl analysis.
	# gene_locs = gene_locs[,c(7,3,4,5)]

	return(gene_locs[,1:4])
} # end get_gene_locs()

#' Takes an ExpressionData object and filters it by the given MAD threshold.
#' @param exp an ExpressionData object
#' @param thres MAD threshold above which to keep genes
#' @return a data.frame with gene locations 
#' @export
mad_filter = function(exp, thres = .05) {
	mad_keep = apply(exp$values, 1, mad) > thres
	exp$values = exp$values[mad_keep,]

	return(exp)
}

#' Takes the path to a DNA methylation matrix and return a MethylationData object.
#' @param meth_data_loc path to a matrix of DNA methylation beta values
#' @return a MethylationData object
#' @export
load_meth = function(meth_data_loc) {
	# Get the methylation data.
	meth_data = MethylationData$new()
	suppressWarnings(meth_data$load_data(meth_data_loc))

	# Handle case of less column names than columns
	num_col = ncol(meth_data$beta_values)
	meth_samples = unlist(strsplit(readLines(meth_data_loc, n=1), "\t"))
	if(num_col == length(meth_samples)) {
		colnames(meth_data$beta_values) = meth_samples
	}
	
	return(meth_data)
}

#' Takes the path to a gene expression matrix and return an ExpressionData object.
#' @param exp_data_loc path to a matrix of gene expression values
#' @return a ExpressionData object
#' @export
load_exp = function(exp_data_loc) {
	# Load expression data
	exp = ExpressionData$new()
	suppressWarnings(exp$load_data(exp_data_loc))

	# Handle case of less column names than columns
	num_col = ncol(exp$values)
	exp_samples = unlist(strsplit(readLines(exp_data_loc, n=1), "\t"))
	if(num_col == length(exp_samples)) {
		colnames(exp$values) = exp_samples
	}

	return(exp)
}