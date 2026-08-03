#' Calculates prognostic statistics based on m2eQTLs.
#'
#' @param meth a MethylationData object, or a data.frame with DNA methylation Beta-values
#' @param exp an ExpressionData object, or a data.frame with gene expression values
#' @param clin a data.frame of clinical data for both meth and exp samples
#' @param eqtls object returned from get_m2eqtls()
#' @param gene_type type of gene list being provided ("trans", "cis", "random", "pre") where trans and cis are those types of m2eQTLs and random is a random set of genes
#' @param gene_sig file to store the gene list in (or load it from)
#' @param covariates a vector of covariates to include from the clinical data
#' @param gene_seed random seed to generate random gene list from
#' @param num_genes number of genes to use for a random gene list
#' @param num_probes number of probes to use for a random probe list
#' @param intersect_data whether to intersect samples for meth and exp (automatically done if integrate_data=TRUE)
#' @param integrate_data whether or not to integrate meth and exp values into a single data frame for the model
#' @param single_data which single data type to use if integrate_data=FALSE, can be "exp" or "meth"
#' @param time_to_event_col name of column with time to event
#' @param event_col name of column with sample event status, or outcome for binomial models
#' @param monte_carlo_size how many random splits of the training data to perform
#' @param training_prop proportion of the data to use for training
#' @param val a list of validation data sets (list of data.frame objects, with combined methylation and expression values or just one or the other)
#' @param val_clin a list of validation clinical data, parallel to val
#' @param alpha a value in [0,1] where 0 is Ridge penalty and 1 is LASSO
#' @return a list of results
#' @export
evaluate = function(
		meth = NULL, 
		exp = NULL, 
		clin = NULL,
		eqtls = NULL, 
		gene_type = NULL,  
		covariates = c(),
		gene_sig = NULL,
		gene_seed = NULL, 
		genes = NULL,
		num_genes = NULL,
		num_probes = NULL,
		intersect_data = FALSE,
		integrate_data = FALSE,
		single_data = "exp",
		time_to_event_col = NULL,
		event_col = NULL,
		rescale = TRUE,
		monte_carlo_size = 100,
		training_prop = .7,
		pam50risk = NULL,
		val = list(),
		val_clin = list(),
		val_pam50 = list(),
		alpha=0) {

	covar_penalty = FALSE

	if(!is.null(time_to_event_col)) {
		type = "cox"
	} else {
		type = "binomial"
	}

	# Set some flags for whether or not to build models for a user 
	# supplied risk score, and whether or not to use independent
	# validation data.
	PAM50 = !is.null(pam50risk)

	# Get gene list to build model from.
	# "trans" or "cis" will get those type of m2eGenes
	# "custom" will simply get a list of genes from a file
	# "random" will return a random gene set
	if(gene_type == "trans" || gene_type == "cis") {
		genes = get_genes(eqtls, gene_type, integrate_data)

		if(!is.null(gene_sig)) {
			write.table(
				genes, gene_sig, sep="\t", col.names=F, row.names=F, quote=F)
		}

		if("MethylationData" %in% class(meth)) {
			meth = meth$m_values[rownames(meth$m_values) %in% get_probes(eqtls),]
		} else {
			meth = meth[rownames(meth) %in% get_probes(eqtls),]
		}
	
	} else if(gene_type == "pre") {
		genes = read.table(gene_sig, sep="\t", header=F, row.names=NULL)
		genes = genes[,1]

		if(integrate_data | intersect_data) {
			sig = genes
			genes = sig[1:num_genes]
			if("MethylationData" %in% class(meth)) {
				meth = meth$m_values[rownames(meth$m_values) %in% sig,]
			} else {
				meth = meth[rownames(meth) %in% sig,]
			}
		}
	} else { # random genes and probes
		if("ExpressionData" %in% class(exp)) {
			genes = get_random_genes(exp$values, gene_seed, num_genes)
		} else {
			genes = get_random_genes(exp, gene_seed, num_genes)
		}
		
		if(integrate_data) {
			if("MethylationData" %in% class(meth)) {
				probes = get_random_genes(meth$m_values, gene_seed, num_probes)
				meth = meth$m_values[rownames(meth$m_values) %in% probes,]
			} else {
				probes = get_random_genes(meth, gene_seed, num_probes)
				meth = meth[rownames(meth) %in% probes,]
			}
		}
	}

	if(!is.null(exp)) {
		# Filter genes with expression data to those identified in m2eQTL analysis
		if("ExpressionData" %in% class(exp)) {
			exp = exp$values[rownames(exp$values) %in% genes,]
		} else {
			exp = exp[rownames(exp) %in% genes,]
		}
	}

	if(rescale) {
		if(!is.null(exp)) {
			# Rescale the exp data from 0 to 1 (by gene).
			scale_info = scale_values(exp)
			exp = scale_info$profile
		}
		if(!is.null(meth)) {
			# Rescale the meth data from 0 to 1 (by probe).
			meth = scale_values(meth)$profile
		}
	}

	if(length(val) > 0) {
		for(i in 1:length(val)) {
			val[[i]] = val[[i]][as.character(rownames(exp)),]

			if(rescale) {
				# Rescale the external validation data.
				val[[i]] = scale_samples(val[[i]], scale_info$means, scale_info$ranges)
				#val[[i]] = scale_values(val[[i]])$profile
			}
		} # end for each validation dataset
	}

	# Transpose the data.
	if(!is.null(exp)) {
		expt = t(exp); rm(exp)
	} else {
		metht = NULL
	}
	if(!is.null(meth)) {
		metht = t(meth); rm(meth)
	} else {
		metht = NULL
	}

	# If the intersect_data flag is set, then subset to samples 
	# with both expression and methylation data.
	if(intersect_data) {
		# Make sure all samples have both expression and methylation data.
		expt = expt[rownames(expt) %in% rownames(metht),]
		metht = metht[rownames(metht) %in% rownames(expt),]
	} 

	# If the integrate_data flag is set, then combine expression
	# values and probes into a single data.frame, else use only
	# one or the other.	
	if(integrate_data) {
		comb = base::merge(expt, metht, by="row.names")
		rownames(comb) = comb[,1]
		comb = comb[,-1]
	} else {
		if(single_data == "exp") {
			comb = expt
		} else {
			comb = metht
		}		
	}

	if(!is.null(expt)) {
		rm(expt)
	}
	if(!is.null(metht)) {
		rm(metht)
	}

	# Subset to clinical data that has molecular data.
	clin = clin[rownames(clin) %in% rownames(comb),]

	if(type == "cox") {
		# Make a survival object for the training data.
		surv = Surv(clin[,time_to_event_col], clin[,event_col])

		# Remove samples without survival data.
		has_surv = !is.na(surv)
		clin = clin[has_surv,]

		# Display number of samples without survival data.
		mm(paste0(sum(!has_surv), 
			" samples dropped for missing survival data."))
	}
	
	# Put molecular data in same order and subset to same samples as clinical 
	# data.
	comb = comb[as.character(rownames(clin)),]

	# Define some lists to collect results.
	pred_list = covar_pred_list = comb_pred_list = list()
	pam50_pred_list = pam50_comb_pred_list = list()
	iv_pred_list = iv_comb_pred_list = list()
	iv_covar_pred_list = iv_pam50_pred_list = list()
	iv_pam50_comb_pred_list = iv2_pred_list = list()
	iv2_pam50_pred_list = iv2_comb_pred_list = list()
	iv2_covar_pred_list = iv2_pam50_comb_pred_list = list()
	ci_list = train_list = test_list = list_of_preds = list()
	list_of_predcombs = list_of_predcovars = list()
	list_of_predpam50s = list_of_predpam50combs = list()
	model_list = nph_list = list()
	lambda_list = list()
	genes = df = c()

	# Variable to hold a running average.
	avg = 0

	# Calculate C-index for N random splits of the data (begin Monte-Carlo)
	# Note, we should probably deal with the missing data before this loop.
	for(i in 1:monte_carlo_size) {

		# For each seed, split data training/test.
		set.seed(i)
		mm(paste0("Seed ", i))
		train = sample(1:nrow(comb), training_prop*nrow(comb), replace=FALSE) 
		train_data = comb[train,]
		test_data = comb[-train,]

		## Subset data to those with molecular and clinical data.
		# This shouldn't be done here, but too late to change it 
		# for now. This also reorders data so samples are in same 
		# order in both clinical and molecular data.
		clin_train = clin[rownames(clin) %in% rownames(train_data),]
		train_data = train_data[as.character(rownames(clin_train)),]

		## Subset data to those with molecular and clinical data.
		clin_test = clin[rownames(clin) %in% rownames(test_data),]
		test_data = test_data[as.character(rownames(clin_test)),]

		# Store the training and test sample ids.
		train_list[[i]] = rownames(clin_train)
		test_list[[i]] = rownames(clin_test)

		if(type == "cox") {
			# Create survival objects for training and test data.
			outcome_train = Surv(
				clin_train[,time_to_event_col], 
				clin_train[,event_col])
			outcome_test = Surv(
				clin_test[,time_to_event_col], 
				clin_test[,event_col])
		} else {
			outcome_train = clin_train[,event_col]
			outcome_test = clin_test[,event_col]
		}

		# Put training data in data.frame.
		train_data = data.frame(train_data)

		# Build M2EFM model.
		#stdz = ifelse(rescale, FALSE, TRUE)
		m2efm = build_m2efm(train_data, outcome_train, clin_train, covariates, 
			standardize=FALSE, alpha=alpha)

		# Put test data in data.frame.	
		test_data = data.frame(test_data)

		if(alpha >= 0) {
			# Get hazard scores for training data.
			pred_train = predict(m2efm$initial_model, 
				data.matrix(train_data),s="lambda.min", type="response")
		} else {
			if(type == "cox") {
				pred_train = predict(m2efm$initial_model, 
					train_data)$predicted
			} else {
				pred_train = predict(m2efm$initial_model, train_data, type="prob")[,1]
			}
		}

		# Subset the user supplied risk scores, if present, to those in 
		# the training data.
		if(PAM50) {
			#pam50_train = pam50risk[as.character(rownames(pred_train)),,
			pam50_train = pam50risk[as.character(rownames(train_data)),,
				drop=FALSE]
		}

		if(length(covariates) > 0) {
			# Get clinical data for the desired covariates for the training data.
			covar_train = clin_train[,c(covariates),drop=F]

			if(type == "cox") {
				# Build covariate-only model.
				if(covar_penalty) {
					covar_train_dat = data.frame(covar_train)
					covar_train_dat = model.matrix(~ ., data=covar_train_dat)
					covar_fit = combfit = cv.glmnet(covar_train_dat, 
						outcome_train, 
						family=type, 
						nfolds=10, 
						standardize=TRUE, 
						alpha=0)
				} else {
					covar_fit = coxph(outcome_train ~ ., data=covar_train)
				}

			} else {
				covar_fit = glm(outcome_train ~ ., data=covar_train, family="binomial")
			}

		}

		val_pred = list()
		val_test = list()

		# Predict the ridge regression model on the validation data.
		if(length(val) > 0) {
			for(j in 1:length(val)) {
				if(alpha >= 0) {
					val_pred[[j]] = predict(m2efm$initial_model, 
						data.matrix(t(val[[j]])),
						s="lambda.min", type="response")
				} else {
					rownames(val[[j]]) = chartr("-", ".", rownames(val[[j]]))
					if(type=="cox") {
						val_pred[[j]] = predict(m2efm$initial_model,
							data.frame(t(val[[j]])))$predicted
					} else {
						rownames(val[[j]]) = chartr("-", ".", rownames(val[[j]]))
						val_pred[[j]] = predict(m2efm$initial_model, 
							t(val[[j]]), type="prob")[,1]
					}

				}
				
				if(length(covariates) > 0) {
					val_test[[j]] = data.frame(pred=val_pred[[j]],
						val_clin[[j]][colnames(val[[j]]),covariates,
						drop=FALSE])
					colnames(val_test[[j]])[1] = "pred"
				}
			}
		}

		# Create a model for PAM50.
		if(PAM50 && length(covariates) > 0) {
			pam50_comb_train = cbind(pred=pam50_train[,1], covar_train)

			colnames(pam50_comb_train) = c("pred", covariates)

			vars = paste0(c("pred", covariates), collapse="+")
			form = as.formula(paste0("outcome_train ~ ", vars))

			if(type == "cox") {
				pam50_comb_fit = coxph(form, data=pam50_comb_train)
			} else {
				pam50_comb_fit = glm(form, data=pam50_comb_train, family="binomial")
			}
			
		}

		if(alpha >= 0) {
			# Predict the ridge regression model on the test data.
			pred = predict(m2efm$initial_model, 
				data.matrix(test_data[,colnames(train_data)]), s="lambda.min", type="response")
		} else {
			if(type=="cox") {
				pred = predict(m2efm$initial_model,
					test_data[,colnames(train_data)])$predicted
			} else {
				pred = predict(m2efm$initial_model, 
					test_data[,colnames(train_data)], type="prob")[,1]
			}
			
		}

		list_of_preds[[i]] = pred

		if(length(covariates) > 0) {
			# Create a combined dataset with the test data risk scores and 
			# covariates.
			covar_test = clin_test[,c(covariates),drop=F]
			comb_test = cbind(pred, covar_test)
			colnames(comb_test) = c("pred", covariates)

			# Predict the final regression model on the test data 
			#(both combined and covariate only).
			if(covar_penalty) {
				comb_test = data.frame(comb_test)
				comb_test = model.matrix(~ ., data=comb_test)
				comb_pred = predict(m2efm$final_model, comb_test, s="lambda.min", type="response")
				covar_test_dat = data.frame(covar_test)
				covar_test_dat = model.matrix(~.,data=covar_test_dat)
				covar_pred = predict(covar_fit, covar_test_dat, s="lambda.min", type="response")
			} else {
				comb_pred = predict(m2efm$final_model, comb_test)
				covar_pred = predict(covar_fit, covar_test)
			}


			list_of_predcovars[[i]] = covar_pred
			list_of_predcombs[[i]] = comb_pred
		}
		
		if(PAM50) {
			# We can get the concordance index from the PAM50 ROR directly.
			# pam50_test = pam50risk[as.character(rownames(pred)),,
			pam50_test = pam50risk[as.character(rownames(test_data)),,
				drop=FALSE]

			if(type == "cox") {
				pam50_pred_list[[i]] = concordance.index(pam50_test[,1], 
					outcome_test[,1], outcome_test[,2])
			} else {
				pam50_roc = AUC::roc(pam50_test[,1], outcome_test)
				pam50_pred_list[[i]] = auc(pam50_roc)
			}

			list_of_predpam50s[[i]] = pam50_test[,1,drop=FALSE]

			if(length(covariates) > 0) {
				# We also combine PAM50 with covariates for an additional
				# model.
				pam50_comb_test = cbind(pam50_test[,1], covar_test)
				colnames(pam50_comb_test) = c("pred", covariates)
				pam50_comb_pred = predict(pam50_comb_fit, pam50_comb_test)

				if(type == "cox") {
					pam50_comb_pred_list[[i]] = concordance.index(
						pam50_comb_pred, outcome_test[,1], outcome_test[,2])
				} else {
					pam50_comb_roc = AUC::roc(pam50_comb_pred, outcome_test)
					pam50_comb_pred_list[[i]] = auc(pam50_comb_roc)
				}

				list_of_predpam50combs[[i]] = pam50_comb_pred
			}
		}

		val_covar_pred = list()
		val_comb_pred = list()

		if(i == 1) {
			if(length(val) > 0) {
				for(j in 1:length(val)) {
					iv_pred_list[[j]] = list()
					iv_covar_pred_list[[j]] = list()
					iv_comb_pred_list[[j]] = list()
					iv_pam50_comb_pred_list[[j]] = list()
					iv_pam50_pred_list[[j]] = list()
				}
			}
		}

		if(length(val) > 0) {
			for(j in 1:length(val)) {
				if(length(covariates) > 0) {
					val_covar_pred[[j]] = predict(covar_fit, val_test[[j]][,covariates,drop=FALSE])
					val_comb_pred[[j]] = predict(m2efm$final_model, val_test[[j]])

					if(type == "cox") {
						iv_covar_pred_list[[j]][[i]] = concordance.index(val_covar_pred[[j]],
							val_clin[[j]][as.character(colnames(val[[j]])), time_to_event_col],
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						iv_comb_pred_list[[j]][[i]] = concordance.index(val_comb_pred[[j]],
							val_clin[[j]][as.character(colnames(val[[j]])), time_to_event_col],
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
					} else {
						iv_covar_roc = AUC::roc(val_covar_pred[[j]],
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						iv_covar_pred_list[[j]][[i]] = auc(iv_covar_roc)

						iv_comb_roc = AUC::roc(val_comb_pred[[j]],
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						iv_comb_pred_list[[j]][[i]] = auc(iv_comb_roc)
					}
					
					
				}

				if(type == "cox") {
					iv_pred_list[[j]][[i]] = concordance.index(val_pred[[j]],
						val_clin[[j]][as.character(colnames(val[[j]])), time_to_event_col],
						val_clin[[j]][as.character(colnames(val[[j]])), event_col])
				} else {
					iv_roc = AUC::roc(val_pred[[j]],
						val_clin[[j]][as.character(colnames(val[[j]])), event_col])
					iv_pred_list[[j]][[i]] = auc(iv_roc)
				}
				

				if(length(val_pam50) > 0) {
					val_pam50_test = cbind(pred=val_pam50[[j]][as.character(rownames(val_test[[j]])),1],
						val_test[[j]][,covariates,drop=FALSE])

					if(type == "cox") {
						iv_pam50_pred_list[[j]][[i]] = concordance.index(val_pam50_test$pred,
							val_clin[[j]][as.character(colnames(val[[j]])), time_to_event_col],
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
					} else {
						val_pam50_roc = AUC::roc(val_pam50_test$pred,
							val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						iv_pam50_pred_list[[j]][[i]] = auc(val_pam50_roc)
					}
					
					if(length(covariates) > 0) {
						val_pam50_comb_pred = predict(pam50_comb_fit, val_pam50_test)

						if(type == "cox") {
							iv_pam50_comb_pred_list[[j]][[i]] = concordance.index(val_pam50_comb_pred,
								val_clin[[j]][as.character(colnames(val[[j]])), time_to_event_col],
								val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						} else {
							val_pam50_comb_roc = AUC::roc(val_pam50_comb_pred,
								val_clin[[j]][as.character(colnames(val[[j]])), event_col])
						}
						
					}
				}
			}	
		} # end for each validation dataset

		# Collect c-indices for the various data.
		if(length(covariates) > 0) {
			if(type == "cox") {
				comb_pred_list[[i]] = concordance.index(
				comb_pred,
				outcome_test[,1], 
				outcome_test[,2])
				covar_pred_list[[i]] = concordance.index(
				covar_pred,
				outcome_test[,1], 
				outcome_test[,2])	
			} else {
				comb_pred_roc = AUC::roc(comb_pred, outcome_test)
				comb_pred_list[[i]] = auc(comb_pred_roc)
				covar_pred_roc = AUC::roc(covar_pred, outcome_test)
				covar_pred_list[[i]] = auc(covar_pred_roc)
			}

		}

		if(type == "cox") {
			pred_list[[i]] = concordance.index(
				pred,
				outcome_test[,1], 
				outcome_test[,2])
		} else {
			pred_roc = AUC::roc(pred, outcome_test)
			pred_list[[i]] = 1 - auc(pred_roc)
		}
		
		model_list[[i]] = m2efm$final_model
		lambda_list[[i]] = m2efm$initial_model$lambda.min
		
		# Display the results for each split of the data.
		if(length(covariates) > 0){
			if(type == "cox") {
				mm(paste(
				comb_pred_list[[i]]$c.index, 
				covar_pred_list[[i]]$c.index, 
				pred_list[[i]]$c.index))	
			} else {
				mm(paste(
					comb_pred_list[[i]],
					covar_pred_list[[i]],
					pred_list[[i]]))
			}
			
		} else {
			if(type == "cox") {
				mm(paste(pred_list[[i]]$c.index))
			} else {
				mm(paste(pred_list[[i]]))
			}
		}

		if(type == "cox" && !covar_penalty) {
			# Keep track of possible violations of proportional hazards assumption.
			zph = cox.zph(m2efm$final_model, transform="rank")
			nph_list[[i]] = rownames(zph$table[,3,drop=FALSE])[zph$table[,3] < .05]
		}

		if(length(covariates) > 0) {
			if(type == "cox") {
				avg = (avg * (i-1) + comb_pred_list[[i]]$c.index) / i
			} else {
				avg = (avg * (i-1) + comb_pred_list[[i]]) / i
			}
			
			mm(paste0("Mean C-index: ", avg))
		}
		
	} # end Monte-Carlo

	# The rest of the code in this function pulls out the C-index 
	# from the concordance results for various combinations of the data.
	pred_res = sapply(pred_list, `[[`, 1)
	
	iv_res_list = list()

	if(length(covariates) > 0) {
		comb_res = sapply(comb_pred_list, `[[`, 1)
		covar_res = sapply(covar_pred_list, `[[`, 1)

		iv_comb_res_list = list()
		iv_covar_res_list = list()
	}

	if(length(val) > 0) {
		for(i in 1:length(val)) {
			iv_res_list[[i]] = sapply(iv_pred_list[[i]], `[[`, 1)
		}
	}
	if(length(val) > 0 && length(covariates) > 0) {
		for(i in 1:length(val)) {
			iv_comb_res_list[[i]] = sapply(iv_comb_pred_list[[i]], `[[`, 1)
			iv_covar_res_list[[i]] = sapply(iv_covar_pred_list[[i]], `[[`, 1)
		}
	}

	evaluate_results = list(
			pred_res = pred_res,
			lambdas = lambda_list)
	if(length(val) > 0) {
		evaluate_results = c(evaluate_results, list(
			iv_pred_res=iv_res_list))
	}

	if(length(covariates) > 0){
		evaluate_results = c(evaluate_results, list(
			comb_res=comb_res,
			covar_res=covar_res))
		if(length(val) > 0) {
			evaluate_results = c(evaluate_results, list(
				iv_comb_res=iv_comb_res_list, 
				iv_covar_res=iv_covar_res_list))
		} 
	}

	if (gene_type %in% c("trans", "cis", "nca")){
		evaluate_results = c(evaluate_results, list(
			train_list=train_list,
			test_list=test_list,
			list_of_preds=list_of_preds,
			comb=comb,
			clin=clin,
			models=model_list,
			nph=nph_list))
		if(length(covariates) > 0) {
			evaluate_results = c(evaluate_results, list(
				list_of_predcombs=list_of_predcombs,
				list_of_predcovars=list_of_predcovars))
		}	
	}

	if(PAM50 && gene_type!="random") {
		pam50_res = sapply(pam50_pred_list, `[[`, 1)
		
		evaluate_results = c(evaluate_results, list(
			pam50_res=pam50_res,
			list_of_predpam50s=list_of_predpam50s))

		if(length(val) > 0) {
			iv_pam50_res_list = list()
			for(i in 1:length(val)) {
				iv_pam50_res_list[[i]] = sapply(iv_pam50_pred_list[[i]], `[[`, 1)
			}

			evaluate_results = c(evaluate_results, list(
				iv_pam50_res=iv_pam50_res_list))
		}	
	} 

	if(PAM50 && gene_type!="random" && length(covariates) > 0) {
		pam50_comb_res = sapply(pam50_comb_pred_list, `[[`, 1)
		evaluate_results = c(evaluate_results, list(
			pam50_comb_res=pam50_comb_res, 
			list_of_predpam50combs=list_of_predpam50combs))

		if(length(val) > 0) {
			iv_pam50_comb_res_list = list()
			for(i in 1:length(val)) {
				iv_pam50_comb_res_list[[i]] = sapply(iv_pam50_comb_pred_list[[i]], `[[`, 1)
			}
			evaluate_results = c(evaluate_results, list(
				iv_pam50_comb_res=iv_pam50_comb_res_list))
		}
	} 

	class(evaluate_results) = "M2EFM_evaluation"
	return(evaluate_results)

} # end evaluate()

#' Build the M2EFM model.
#'
#' @param data data.frame with molecular data. 
#' @param outcome Survival object from survival package, or binomial outcome
#' @param clin data.frame with clinical data
#' @param covariates vector of column names in clinical data
#' @param standardize whether or not to standardize the data
#' @return a vector of genes to use for models
#' @export
build_m2efm = function(data, outcome, clin, covariates, standardize=FALSE, alpha=0) {
	type = ifelse("Surv" %in% class(outcome), "cox", "binomial")

	# Whether or not to penalize covariates.
	covar_penalty = FALSE

	if(alpha >= 0) {
		# Build initial model
		glmfit = cv.glmnet(data.matrix(data), 
			outcome, 
			family=type, 
			nfolds=10, 
			standardize=standardize, 
			alpha=alpha)
				
		# Get hazard scores for data.
		pred_data = predict(glmfit, data.matrix(data),s="lambda.min", type="response")
	} else {
		library(randomForestSRC)

		#t_data = data.frame(t(data))

		if(type == "cox") {
			data = cbind(data, time = outcome[,1], status=outcome[,2])

			# Build initial model
			#rffit = rfsrc(Surv(time, status)~., data=data, ntree=1000, na.action="na.impute", splitrule="logrank", nsplit=1, importance="random", seed=-1)
			if(alpha == -2) {
				rffit = rfsrc(Surv(time, status)~., data=data, ntree=1000, na.action="na.impute", splitrule="logrank", nsplit=1, importance="random", seed=-1)
			} else {
				rffit = rfsrc(Surv(time, status)~., data=data, ntree=1000, na.action="na.impute", nsplit=1, importance="random", seed=-1)
			}

			# Get hazard scores for data.
			pred_data = predict(rffit, data)$predicted
			
		} else {
			rffit = randomForest(data, outcome, ntree=1000)

			pred_data = predict(rffit, data, type="prob")[,1]
		}
		
		glmfit = rffit
	}

	if(length(covariates) > 0) {
		# Get clinical data for the desired covariates for the data.
		covar = clin[,c(covariates),drop=F]
			
		# Create a combined dataset with the training data risk scores and the 
		# covariates.
		comb = cbind(pred=pred_data, covar)
		colnames(comb) = c("pred", covariates)

		vars = paste0(c("pred", covariates), collapse="+")
	} else {
		# If no covariates, build the model without them.
		comb = data.frame(pred=pred_data)
		colnames(comb)[1] = "pred"

		vars = "pred"
	}

	form = as.formula(paste0("outcome ~ ", vars))

	# Build final model.
	if(type == "cox") {
		if(covar_penalty) {
			comb_dat = model.matrix(outcome ~ ., data=comb)
			combfit = cv.glmnet(comb_dat, 
			outcome, 
			family=type, 
			nfolds=10, 
			standardize=standardize, 
			alpha=alpha)
		} else {
			combfit = coxph(form, data=comb)
		}

	} else {
		combfit = glm(form, data=comb, family="binomial")
	}
	
	m2efm_result = list(initial_model=glmfit, final_model=combfit)
	class(m2efm_result) = "M2EFM_models"

	return(m2efm_result)
} # end build_m2efm

#' Returns the genes to use from m2eQTL analysis.
#'
#' @param eqtls object returns from get_m2eqtls()
#' @param gene_type "cis" or "trans"
#' @return a vector of genes to use for models
#' @export
get_genes = function(eqtls, gene_type, integrate_data=FALSE) {
	if(gene_type=="trans") {
		genes = unique(eqtls$trans_link$gene)
		more_genes = unique(eqtls$result$cis$eqtls$gene[eqtls$result$cis$eqtls$snps %in% eqtls$result$trans$eqtls$snps])
		#more_genes = unique(eqtls$result$cis$eqtls$gene)
		if(!integrate_data) {
			genes = unique(c(as.character(genes), as.character(more_genes)))
		}
	} else {
		genes = unique(eqtls$cis_link$gene)
	}

	return(genes)
} # end get_genes()

#' Returns the probes to use from m2eQTL analysis.
#'
#' @param eqtls object returns from get_m2eqtls()
#' @return a vector of probes to use for models
#' @export
get_probes = function(eqtls) {
	return(as.character(unique(eqtls$trans_link$snps)))
} # end get_probes()

#' Returns the genes to use selected at random.
#'
#' @param exp data.frame of gene expression data
#' @param seed a random seed
#' @param num_genes the number of genes to take at random
#' @return a vector of genes to use for models
#' @export
get_random_genes = function(exp, seed, num_genes) {
	set.seed(seed)
	genes = sample(1:nrow(exp), num_genes, replace=FALSE)
	genes = rownames(exp)[genes]

	return(genes)
} # end get_random_genes()

#' Converts AJCC stages to a simple four stage form.
#'
#' @param x a vector of AJCC stages (e.g. Stage IA)
#' @return a vector of simplified AJCC stages
#' @export
reduce_stage = function(x) {
	'Reduces a string containg cancer stage to simple 4 stage factor.
	'
	# Collapse tumor stages.
	x[x %in% c("Stage I", "Stage IA", "Stage IB", "I", "IA", "IB")] = "Stage I"
	x[x %in% c("Stage II", "Stage IIA", "Stage IIB", "II", "IIA", "IIB")] = 
		"Stage II"
	x[x %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC", "III", 
		"IIIA", "IIIa", "IIIB", "IIIC")] = "Stage III"

	if(is.factor(x)) {
		# Get rid of empty levels so they aren't used in regression.
		x = droplevels(x)
	}

	return(x)
} # reduce_stage()

#' Get pam50 rorS scores for sample with gene expression data. It
#' should go without saying that these should be breast cancer samples.
#'
#' @param exp a data.frame with gene expression values
#' @return a list with subtype calls and rorS scores for the samples in exp
#' @export
get_pam50 = function(exp) {
	load_it(c("org.Hs.eg.db", "genefu"))

	x = org.Hs.eg.db::org.Hs.egALIAS2EG
	mapped_genes = mappedkeys(x)
	xx = AnnotationDbi::as.list(x[mapped_genes])

	exp_xx = exp[rownames(exp) %in% names(xx),]
	entrez_ids = sapply(xx[rownames(exp_xx)],function(x) x[1])

	annot = data.frame(GeneSymbol = rownames(exp_xx), 
		EntrezGene.ID = entrez_ids)

	exp_t_xx = t(exp_xx)

	RORs = rorS(data=exp_t_xx, annot=annot)
	pam50risk = as.data.frame(RORs$score)
	subtypes = molecular.subtyping(sbt.model = "pam50",
		data = exp_t_xx,annot=annot) 
	subtype_calls = subtypes$subtype
	names(subtype_calls) = chartr(".", "-", names(subtype_calls))

	rownames(pam50risk) = chartr(".", "-", rownames(pam50risk))

	return(list(subtype_calls=subtype_calls, pam50risk=pam50risk))
} # end get_pam50()

#' Center and scale the expression values.
#'
#' @param profile a data.frame with gene expression values
#' @return a data.frame of gene expression values with each gene in [-1,1]
#' @export
scale_values = function(profile) {
	rowmins = matrixStats::rowMins(data.matrix(profile))
	rowmaxes = matrixStats::rowMaxs(data.matrix(profile))
	rowmeans = rowMeans(profile)
	ranges = rowmaxes - rowmins

	profile = scale_samples(profile, rowmeans, ranges)

	return(list(profile=profile, means=rowmeans, ranges=ranges))
}

#' Scales individual samples to another dataset.
#'
#' @param profile a data.frame with gene expression values
#' @param means a vector of means by gene from original dataset
#' @param ranges a vector of ranges or sds by gene from original dataset
#' @return a data.frame of gene expression values with each gene scaled
#' @export
scale_samples = function(profile, means, ranges) {
	for(i in 1:ncol(profile)) {
		profile[,i] = (profile[,i] - means)/ranges
	}
	return(profile)
}

#' Display a message to the user.
#'
#' @param text text of the message to display
#' @return nothing
#' @export
mm = function(text, appendLF=TRUE) {
	message(text, appendLF=appendLF)
}

#' Summarize the results of an M2EFM models.
#'
#' @param m2efm_result object returned by build_m2efm()
#' @return nothing
#' @export
summary.M2EFM_models = function(m2efm_result, ...) {
	summary(m2efm_result$final_model)
} # end summary.M2EFM_models


#' Create boxplots of the evaluation results.
#'
#' @param eval_result object returned by evaluate()
#' @return nothing
#' @export
plot.M2EFM_evaluation = function(eval_result, pal=NULL, ...) {
	load_it("ggplot2")

	result = data.frame(`C-index`=c(
		eval_result$covar_res,
		eval_result$comb_res,
		eval_result$pred_res),
		`Model`=factor(rep(c(
			"Covariates",
			"M2EFM+Covariates", 
			"M2EFM"), each=100), 
		levels=c(
		"Covariates",
		"M2EFM+Covariates",
		"M2EFM")))
	
	if("iv_comb_res" %in% attributes(eval_result)$names &&
		length(eval_result$iv_comb_res) > 0) {
		for(i in 1:length(eval_result$iv_comb_res)) {
			row = data.frame(`C-index`=c(
				eval_result$iv_covar_res[[i]],
				eval_result$iv_comb_res[[i]],
				eval_result$iv_pred_res[[i]]),
				`Model`=factor(rep(c(
					paste0("Validation ", i, " Covariates"),
					paste0("Validation ", i, " M2EFM+Covariates"), 
					paste0("Validation ", i, " M2EFM")), each=100),
				levels=c(
					paste0("Validation ", i, " Covariates"),
					paste0("Validation ", i, " M2EFM+Covariates"), 
					paste0("Validation ", i, " M2EFM"))))
			result = rbind(result, row)
		}
	}

	if(!is.null(pal)) {
		cbPalette = pal
	} else {
		cbPalette = rep("#FFFFFF", times=length(levels(result$Model)))
	}
	
	ggplot(result) + 
		geom_boxplot(aes(x=Model, y=C.index, fill=Model), notch=TRUE) +
		theme_bw() + xlab("Model") + ylab("C-index") +
		scale_fill_manual(values=cbPalette) + 
		theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1, size=14),
			axis.title.x = element_text(size=14, face="bold")) +
		theme(axis.text.y = element_text(size=14),
			axis.title.y = element_text(size=14, face="bold")) +
		# facet_grid(~ Dataset, scales="free") +
		theme(strip.text.x = element_text(size = 10, face="bold"), 
			strip.background = element_rect(color="black", fill="white"))
} # end plot.M2EFM_evaluation

#' Summarize the results of an M2EFM evaluation.
#'
#' @param eval_result object returned by evaluate()
#' @return nothing
#' @export
summary.M2EFM_evaluation = function(eval_result, nph=FALSE, ...) {
	if("comb_res" %in% attributes(eval_result)$names &&
		length(eval_result$comb_res) > 0) {
		median_results = data.frame(`Median Result`=c(median(eval_result$covar_res),
			median(eval_result$comb_res), median(eval_result$pred_res)),
		Method=c("Covariates", "M2EFM+Covariates", "M2EFM"))
	} else {
		median_results = data.frame(`Median Result`=median(eval_result$pred_res),
			Method="M2EFM")
	}

	if("iv_pred_res" %in% attributes(eval_result)$names && 
		length(eval_result$iv_pred_res) > 0) {
		if("iv_comb_res" %in% attributes(eval_result)$names &&
			length(eval_result$iv_comb_res) > 0) {
				for(i in 1:length(eval_result$iv_comb_res)) {
					median_results = rbind(median_results, 
						data.frame(`Median Result`=c(
							median(eval_result$iv_covar_res[[i]]),
							median(eval_result$iv_comb_res[[i]]), 
							median(eval_result$iv_pred_res[[i]])),
						Method=c(paste0("Validation ", i, " Covariates"), 
							paste0("Validation ", i, " M2EFM+Covariates"),
							paste0("Validation ", i, " M2EFM"))))
				}
		} else {
			for(i in 1:length(eval_result$iv_res)) {
				median_results = rbind(median_results, 
					data.frame(`Median C-index`= 
						median(eval_result$iv_pred_res[[i]]),
					Method= paste0("Validation ", i, " M2EFM")))
			}
		}

	}

	med_res = median_results
	rownames(med_res) = med_res$Method
	med_res = med_res[,-2,drop=FALSE]
	print(med_res)

	if(nph) {
		cat(sprintf("\nPossible violations of proportional hazards assumption:\n"))
		nph_viol = plyr::count(unlist(eval_result$nph))
		rownames(nph_viol) = nph_viol$x
		nph_viol = nph_viol[,-1,drop=FALSE]
		print(nph_viol)
	}
} # end summary.M2EFM_evaluation


#' Predict the M2EFM model on new data.
#'
#' @param fit the M2EFM models to predict - the result of build_m2efm
#' @param new_data data.frame with molecular data.
#' @param clin data.frame with clinical data
#' @param covariates vector of column names in clinical data
#' @return a data.frame of up to two columns of predictions: for molecular and covariates and molecular only
#' @export
predict.M2EFM_models = function(fit, data, clin, covariates) {
	prediction = predict(fit$initial_model, data.matrix(data), s="lambda.min", type="response")

	if(length(covariates) > 0) {
		new_data = data.frame(pred=prediction, 
			clin[rownames(data), covariates, drop=FALSE])
		colnames(new_data)[1] = "pred"

		comb_prediction = predict(fit$final_model, new_data)
	}

	return(data.frame(mol_and_covar=comb_prediction, mol_only=prediction))
} # end predict.M2EFM_models







