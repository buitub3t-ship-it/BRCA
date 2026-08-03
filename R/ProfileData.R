#'A class for holding gene or methylation profile data.
#'
#'@name ProfileData
#'@title ProfileData
#'@field values The data values.
#'@export ProfileData
#'@exportClass ProfileData
ProfileData <- R6Class("ProfileData",
	#inherit = ProfileData,
	private = list(
	#		
			),
	public = list(
		values = NULL,
		
		initialize =  function(values) {
			require(data.table)
			require(impute)

			if(hasArg(values)) {
				self$values = values
			}
		}, # end initialize method
		
		load_data = function(filename) {
		  'Load beta or expression values from filename.'

		  values = suppressWarnings(fread(filename, data.table=F))
		  rownames(values) = values[,1]
		  values = values[,-1]

		  self$values = values
		  invisible(values)
		}, # end load_data

		impute = function(na_max_prop=.5, seed=1, k=10) {

			# Filter probes more than half NAs
			self$values = self$values[
				rowMeans(is.na(self$values)) <= na_max_prop,]

			self$values = data.frame(impute.knn(
				as.matrix(self$values), 
				k=k, 
				rng.seed=seed)$data)

			colnames(self$values) = chartr(".", "-", colnames(self$values))

			invisible(self$values)
		} # end impute
		
	) #end public
) #end ProfileData
