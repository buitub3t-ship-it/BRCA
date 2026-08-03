query_genes_by_coordinates = function(coordinates, dist) {
	# Given a set of coordinates and a distance, return a data.frame 
	# of gene information for genes that occur within dist of those 
	# coordinates.
	#
	#>coordinates -- a character string in the form of: 
	# chr<num>:<start>-<end>
	#>dist -- distance in base-pairs
	#<returns a vector containing gene symbols
	load_it(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "RSQLite", "stringr"))

	# Get the coordinates of all CpGs in the data.
	loc_matches = extract_locations(coordinates)

	conn = dbConnect(
		RSQLite::SQLite(), 
		dbname=paste0(find.package(
			"TxDb.Hsapiens.UCSC.hg19.knownGene"),
			"/extdata/TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite"))
	coords = dbGetPreparedQuery(conn, 
		paste0("select distinct gene.gene_id, cds._cds_id, cds_chrom, ",
			"cds_start, cds_end, cds_strand ",
			"from cds, gene, splicing where ",
			"gene._tx_id=splicing._tx_id and ",
			"splicing._cds_id=cds._cds_id and ",
			"cds_chrom = ? and abs(cds_start-?) <= ? and",
			" abs(cds_end-?) <= ?"),
		data.frame(
			chr=loc_matches$chrom, 
			loc=loc_matches$start, 
			dist, 
			loc2=loc_matches$end, 
			dist))
	dbDisconnect(conn)

	symbols = query_gene_symbols_by_id(coords[,1])

	out = data.frame(gene_id=character(), symbol=character())

	for(i in 1:nrow(coords)) {
		if(coords[i,1] %in% symbols[,1]) {
			out = rbind(out, unique(symbols[symbols[,1] == coords[i,1],]))
		} else {
			out = rbind(out, data.frame(gene_id=coords[i,1], symbol=""))
		}
	}

	coords = cbind(out, coords)
	return(coords)
} #end query_genes_by_coordinates

query_gene_symbols_by_id = function(entrez_ids) {
	# Given a vector of Entrez ids, return a data.frame containing gene
	# symbols.
	#>entrez_ids -- a vector containing Entrez ids
	#<returns a data.frame containing gene symbols

	load_it("org.Hs.eg.db")

	conn = dbConnect(RSQLite::SQLite(), 
		dbname=paste0(find.package("org.Hs.eg.db"),
			"/extdata/org.Hs.eg.sqlite"))
	symbols = dbGetPreparedQuery(conn,
		paste0('select distinct genes.gene_id, gene_info.symbol from gene_info, genes', 
				' where genes.gene_id = ? and',
				' genes._id=gene_info._id'),
		data.frame(entrez_ids))
	dbDisconnect(conn)

	return(symbols)
} #end query_gene_symbols_by_id

query_gene_ids_by_symbol = function(symbols) {
	load_it(c("RSQLite", "org.Hs.eg.db"))

	# Find gene ids for genes by their symbols.
	x = org.Hs.eg.db::org.Hs.egALIAS2EG
	mapped_genes = mappedkeys(x)
	xx = as.list(x[mapped_genes])

	found_ids = symbols %in% names(xx)
	entrez_ids = sapply(xx[symbols[found_ids]], function(x) x[1])
	
	return(entrez_ids)
} #end query_gene_id_by_synbol

query_gene_coordinates = function(entrez_ids, min_start, max_end) {
	# Given a vector of Entrez ids, return a data.frame containing gene
	# coordinates.
	#>entrez_ids -- a vector containing Entrez ids
	#>min -- Boolean indicating only min starts should be returned
	#>max -- Boolean indicating only max ends should be returned
	#<returns a data.frame containing gene coordinates
	load_it(c("TxDb.Hsapiens.UCSC.hg19.knownGene", "RSQLite"))

	if(min_start) {
		cds_start = 'min(cds_start)'
	} else {
		cds_start = 'cds_start'
	}
	if(max_end){
		cds_end = 'max(cds_end)'
	} else {
		cds_end = 'cds_end'
	}
	
	conn = dbConnect(RSQLite::SQLite(), 
		dbname=paste0(find.package("TxDb.Hsapiens.UCSC.hg19.knownGene"),
		"/extdata/TxDb.Hsapiens.UCSC.hg19.knownGene.sqlite"))
	
	coords =
			dbGetPreparedQuery(conn, 
				paste0('select distinct gene.gene_id, cds._cds_id, cds_chrom, ',
						cds_start, ', ', 
						cds_end, ', ', 
						'cds_strand from cds, gene, splicing where gene.gene_id=? and ',
						'gene._tx_id=splicing._tx_id and splicing._cds_id=cds._cds_id'),
				data.frame(entrez_ids))
	dbDisconnect(conn)
	
	return(coords)
} #end query_gene_coordinates