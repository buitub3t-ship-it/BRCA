#'A class for holding methylation profile data from Illumina 450k.
#'
#'@name MethylationData450k
#'@title MethylationData450k
#'@export MethylationData450k
#'@exportClass MethylationData450k
MethylationData450k <- R6Class("MethylationData450k",
	inherit = MethylationData,
	public = list(	
		initialize =  function(values) {
			super$initialize(values)
		}, # end initialize method
			
		filter_sex_chrom_probes = function() {
			load_it(c("TxDb.Hsapiens.UCSC.hg19.knownGene",
				"RSQLite",
				"stringr",
				"IlluminaHumanMethylation450k.db"))

			# Get the location of the 450k annotation database.
			db = system.file("extdata", "Illumina450kAnnotation.sqlite", package="M2EFM")

			# Get the name of the table in the 450k annotation db.
			illumina_450k_table = "Ilmn450kv12"

			conn = dbConnect(SQLite(), dbname=db)

			probes = dbGetQuery(conn, 
				paste0("select IlmnID from ",
					illumina_450k_table,
					" where CHR = 'X' or CHR = 'Y' "))
			dbDisconnect(conn)

			return(self$values[!rownames(self$values) %in% probes[,1],])
		}, #end filter_sex_chrom_probes

		filter_probes_with_snps = function(dist=10) {
			load_it("Illumina450ProbeVariants.db")

			e = new.env()
			data(probe.450K.VCs.af, envir=e)

			# filter probes with any variants within 10 or 50 bp of the CpG
			if(dist==10) {
				probes = rownames(
					e$probe.450K.VCs.af[e$probe.450K.VCs.af$probe10VC.AMR > 0,])
				return(self$values[!rownames(self$values) %in% probes,])
			} else {
				probes = rownames(
					e$probe.450K.VCs.af[e$probe.450K.VCs.af$probe50VC.AMR > 0,])
				return(self$values[!rownames(self$values) %in% probes,])
			}
			
		}, #end filter_probes_with_snps

		filter_cross_reactive_probes = function() {
			# Returns a list of probes that are cross-reactive.
			# > contraints -- not currently used.
			# < return nothing

			# Get the location of the cross-reactive probes database.
			db = system.file("extdata", "crossReactive.txt", package="M2EFM")

			# Read the list of probes.
			probes = fread(db, header=F, data.table=F) 
			
			return(self$values[!rownames(self$values) %in% probes[,1],])
		}, #end filter_cross_reactive_probes

		collapse_to_islands = function(groups=NULL) {
			'Collapse the beta values to be the mean values over the islands they are 
			in. 
					
			beta_values A data.frame() containing beta values. The rows are the
			probes and the columns are the samples.
					
			returns A data.frame() containing the beta values for which the rows
			are now CpG islands.'

			if(is.null(groups)) {
				db = system.file("extdata", "CGI.sqlite", package="M2EFM")
				
				# Open a connection to the db.
				conn = dbConnect(SQLite(), dbname=db)
				
				# Query for all the islands.
				cpg_islands = dbGetQuery(conn, "select * from Islands")
				
				dbDisconnect(conn)

				groups = cpg_islands[cpg_islands$Probe_ID %in% 
					rownames(self$values),]

				# Get the mean of each group in the data,
				# (CpG islands).
				beta.cgi = aggregate(self$values[as.character(groups$Probe_ID),],
					by=list(groups$cginame),mean)
			} else {
				# Get the mean of each group in the data,
				beta.cgi = aggregate(self$values,by=list(groups),mean)
			}

			# Set the rownames of the collapsed beta values to be the
			# location of the island.
			rownames(beta.cgi) = beta.cgi[,"Group.1"]

			# Remove the column holding the location of the island and
			# return.
			return(subset(beta.cgi, select=-Group.1))
		} #end collapse_to_islands
	) #end public
) #end MethylationData450k
