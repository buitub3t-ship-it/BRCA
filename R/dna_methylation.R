diff_meth = function(
	data1, 
	data2, 
	cutoff=.05, 
	covariates=NULL, 
	clin=NULL, 
	method="eBayes", 
	lfc=0, 
	ref=2) {
  '
  Perform differential methylation analysis on groups of methylation beta
	values.
			
	data1,data2: data.frames of beta values for which the rows are probes and 
    the columns are samples
  cutoff: cutoff for minimum statistical significance to be returned
  covariates: vector of covariate names from clinical data to use
  clin: if you supply covariates, you must also supply a data.frame with
    clinical data that contain the actual covariates with colnames matching
    the covariates vector
  method: "eBayes" or "treat" for moderating statistics
  lfc: minimum required log-fold change for results
  ref: which data are the reference
    
  returns: object of results
  '
  
  load_it(c("limma", "lumi"))      
  
	prefix1 = "Data1"
	prefix2 = "Data2"

  # Prefix each sample id to show its origin.
  colnames(data1) = paste0(prefix1, colnames(data1))
  colnames(data2) = paste0(prefix2, colnames(data2))

  if(ref == 1) {
  	reference = "Data1"
  	target = "groupsData2"
  } else {
  	reference = "Data2"
  	target = "groupsData1"
  }

  # Create a data.frame with a column named Group, that has the value Data1- or Data2- and
  # rownames of the patient ids from the different groups.
	design = data.frame(
    Group = relevel(
			factor(substr( 
				c(colnames(data1), 
				colnames(data2)),1,5)), 
			ref=reference), 
		row.names=c(colnames(data1), colnames(data2)))

	# Get the groups from the design matrix.
	groups = design$Group

	# If there are covariates, add them to the design. In either case, create a
	# model matrix that includes all sample ids and groups.
	if(!is.null(covariates)) {
		covar = clin[rownames(clin) %in% 
					substr(rownames(design),nchar(prefix1),nchar(colnames(data1)[1])),
              	c(covariates),
              	drop=F]
    form = paste0('~ ', paste0(c(covariates, "groups")))

  	DesMat = model.matrix(as.formula(form), data=cbind(covar, groups))
	
	} else {
		DesMat = model.matrix(~groups)
	}

  # combine all data
  comb_data = cbind(data1, data2)
	
  # Use a linear model to try to separate labelled data.
  DMRfit = lmFit(comb_data, DesMat)

  if(method == "eBayes") {
    # Use empirical Bayes approach to test differential methylation,
    # using the model trained before.
    DMRfitEb = eBayes(DMRfit, robust=TRUE)

    # Get a table of the differentially methylated regions below the cutoff.
    DMR = topTable(DMRfitEb, coef=target, number=Inf, p.value=cutoff)
  } else {
    DMRfitEb = treat(DMRfit, lfc=lfc, robust=TRUE)
    DMR = topTreat(DMRfitEb, coef=target, number=Inf, p.value=cutoff)
  }
  
  return(list(DMRfit=DMRfit,DMRfitEb=DMRfitEb,DMR=DMR))
} #end diff_meth
		
query_450k_annotation = function(columns, constraints=NULL, probe_ids=NULL) {
	# Frequently, the full annotation won't be needed. Just get the
	# columns passed.
	#
	# >columns -- a vector of column names to retrieve
	# >contraints -- such as MAPINFO > 1000000 and CHR = 7
	# >probe_ids -- a vector of probe_ids to keep
	# <returns the contents of the columns asked for
	load_it("stringr")
	load_it("RSQLite")

	# Get the location of the 450k annotation database.
	db = system.file("extdata", "Illumina450kAnnotation.sqlite", package="M2EFM")

	illumina_450k_table = "Ilmn450kv12"

	conn = dbConnect(SQLite(), dbname=db)
	sub_table = NULL

	# Constraints passed by user get built into a where clause.
	if(length(constraints)>0) {
		where_clause = paste0(' where ', constraints)
	} else {
		where_clause = ''
	}

	# Put probe ids passed by user into a data.frame and then
	# into an in-memory database table. Only annotations for
	# these probes will be returned.
	probe_ids = data.frame(IlmnID = probe_ids)
	if(length(probe_ids)>0){
		dbGetQuery(conn, "attach ':memory:' as mem")
		dbWriteTable(conn, "mem.probes", probe_ids)
	}

	# If the user specified probes, select from the annotation
	# based on those probes and, potentially, other user constraints.
	# Otherwise select all annotations that match any user constraints.
	if(length(probe_ids)>0) {
		if(str_length(where_clause) > 0) {
			probe_clause = paste0(" and ", illumina_450k_table, 
					".IlmnID = probes.IlmnID")
		} else {
			probe_clause = paste0(" where ", illumina_450k_table, 
					".IlmnID = probes.IlmnID")
		}
		
		columns = paste0(illumina_450k_table, ".", columns)
		col_str = paste(columns, collapse=", ")
		
		subtable = dbGetQuery(conn, 
				paste0('select ', col_str, 
						' from ', illumina_450k_table,
						', mem.probes ',
						where_clause, probe_clause))
	} else{
		col_str = paste(columns, collapse=", ")
		subtable = dbGetQuery(conn, 
			paste0('select ', col_str, ' from ', illumina_450k_table, where_clause))
	}

	# Disconnet the user probe table.
	if(length(probe_ids)>0) {
		dbGetQuery(conn, 'detach mem')
	}

	# Disconnect the database.
	dbDisconnect(conn)

	# Return the selection from the 450k annotation.
	return(subtable)
} #end query_450k_annotation
