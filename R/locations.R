extract_locations = function(regions) {
	'Takes genomic coordinates as chr<num>:<start>-<end> and returns
	a data.frame with the parts broken out.'
	
	load_it("stringr")

	loc_matches = data.frame(do.call(rbind, str_match_all(regions, 
		"(chr\\d+|chrX|chrY):(\\d+)-(\\d+)")), row.names=1:length(regions))
	loc_matches$X2 = as.character(paste(loc_matches$X2))
	loc_matches$X3 = as.numeric(paste(loc_matches$X3))
	loc_matches$X4 = as.numeric(paste(loc_matches$X4))
	colnames(loc_matches) = c("location", "chrom", "start", "end")

	return(loc_matches)
} # end extract_locations