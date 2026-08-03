#'A class for holding methylation profile data.
#'
#'@name MethylationData
#'@title MethylationData
#'@export MethylationData
#'@exportClass MethylationData
MethylationData <- R6Class("MethylationData",
	inherit = ProfileData,
	private = list(
			value_type = "beta"
			),
	active = list(
	    beta_values = function(value) {
	    	' Create beta_values property. In order to save space and time 
			only one copy of the data is kept and converted between M and beta
			values as needed.
	    	'
			if(missing(value)) {
				if(private$value_type!="beta") {
					self$values = m2beta(self$values)
					private$value_type = "beta"
				}
				return(self$values)
			}
			else self$values = value
	    }, # end beta_values property
 
	    m_values = function(value) {
	    	' Create m_values property. In order to save space and time 
			only one copy of the data is kept and converted between M and beta
			values as needed.
	    	'
	    	if(missing(value)) {
	    		if(private$value_type=="beta") {
	    			self$values = beta2m(self$values)
	    			private$value_type = "mvalue"
	    		}
	    		return(self$values)
	    	} else {
				self$values = value
	    	}
	    } # end m_values property
  	),
	public = list(	
		initialize =  function(values) {
			super$initialize(values)
			require(lumi)
		}, # end initialize method
			
		shift_ucsc = function() {
			'This method simply adds .5 to all beta values in beta_values to reverse 
			the zero-centering done by UCSC.'
			self$values = self$values + .5
			message('Data values shifted + .5')
		} #end shift_ucsc
		
	) #end public
) #end MethylationData
