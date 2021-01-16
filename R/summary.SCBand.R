summary.SCBand <- function(object, ...) 
{
	stopifnot(inherits(object,"SCBmeanfd"))
	callfun <- as.character(object$call[1])
	cat("\nCall:\n", paste(deparse(object$call), sep = "\n", collapse = "\n"), 
		"\n\n", sep = "")

	cat("Data:\n")
	if (is.list(object$x)) {
	df <- data.frame( c(paste0("[", round(min(object$x[[1]]),3),
		", ",round(max(object$x[[1]]),3),"]"),		
		if (!is.null(object$y)) {
			c(paste0("[",round(min(object$y[[1]]),3),", ",
			round(max(object$y[[1]]),3),"]"), 
			nrow(object$y[[1]])) }), 
			c(paste0("[", round(min(object$x[[1]]),3),", ",
			round(max(object$x[[1]]),3),"]"),		
			if (!is.null(object$y)) {
			c(paste0("[",round(min(object$y[[2]]),3),", ",
			round(max(object$y[[2]]),3),"]"), 
			nrow(object$y[[2]]))
		})
	)
	colnames(df) <- paste("Data",1:2)
	rownames(df) <- c("Range.x", if(!is.null(object$y)) c("Range.y", "Sample size"))
	print(df, print.gap = 2L, right = FALSE, digits = 4L, row.names = TRUE)
	}
	else {
	cat("Range(x): [", round(min(object$x),3),", ",round(max(object$x),3),"]\n", sep="")
		if (!is.null(object$y)) {
			cat("Range(y): [",round(min(object$y),3),", ",round(max(object$y),3),"]\n", sep="")
			cat("Sample size:", nrow(object$y),"\n")
		} 
	}
	
	
	cat("\nAnalysis:\n")
	if (callfun == "scb.mean")
	# if (is.null(object$par) && (!is.matrix(object$nonpar)) ) {		
		cat("Mean function estimation\n")
		cat("Bandwidth:", round(object$bandwidth,4), "\n")
		cat("Grid size:", object$gridsize, "\n")
		cat("SCB type:", switch(object$scbtype, normal = "normal", 
			bootstrap = "boostrap", tGKF = "tGKF", 
			all = "normal, bootstrap, tGKF"),"\n")
		# if (object$scbtype != "no") {
		cat("Confidence level:", object$level, "\n")
		cat("Replicates:", switch(object$scbtype, 
			normal = paste(object$nrep,"(normal)\n"), 
					bootstrap = paste(object$nboot,"(bootstrap)\n"), 
					all = paste(object$nrep,"(normal)", 
						object$nboot,"(bootstrap)\n")))
		cat("Quantile used for SCB:\n")		
		statresult <- switch(object$scbtype,
			normal = data.frame(normal = object$qnorm),
			boot = data.frame(bootstrap = object$qboot),
			tGKF = data.frame(tGKF = object$qtGKF),
			all = data.frame(normal = object$qnorm, 
				bootstrap = object$qboot, tGKF = object$qtGKF)
		)
		# statresult <- if (object$scbtype == "normal") { data.frame(object$qnorm)
								# } else if (object$scbtype == "boot") { data.frame(object$qboot)
								# } else if (object$scbtype == "tGKF") { data.frame(object$qtGKF)
								# } else data.frame(object$qnorm, object$qboot, object$qtGKF)
				# names(statresult) <- c(switch(object$scbtype, normal = "normal", 
				# bootstrap = "bootstrap", tGKF = "tGKF", all = c("normal", "bootstrap", "tGKF")))
		print(statresult, print.gap = 2L, right = FALSE, digits = 4L, 
			row.names = FALSE)
		# }

	} 
	else {
		
		pfun <- function(p) {
			if (is.null(p)) return(NA)
			if (p < 1e-16) return("< 1e-16")
			if (p < 1e-4) return(paste0("<1e", 
				as.integer(ceiling(log(p, 10)))))
			return(p) 
		}
		statresult <- switch(object$scbtype,
			normal = data.frame(object$teststat, pfun(object$qnorm)),
			boot = data.frame(object$teststat, pfun(object$qboot)),
			tGKF = data.frame(object$teststat, pfun(object$qtGKF)),
			all = data.frame(object$teststat, pfun(object$qnorm), 
				pfun(object$qboot), pfun(object$qtGKF))
		)
		names(statresult) <- switch(object$scbtype,
			normal = c("stat", "normal p"),
			boot = c("stat", "bootstrap p"),
			tGKF = c("stat", "tGKF p"),
			all = c("stat", "normal p","bootstrap p", "tGKF p")			
		)
		# Pnorm <- if (is.null(object$pnorm)) { 
			# NULL
			# } else if (object$pnorm < 1e-16) { 
			# "<1e-16"
		 	# } else if (object$pnorm < 1e-4) { 
			# paste0("<1e", as.integer(ceiling(log(object$pnorm, 10))))  
		   	# } else object$pnorm

		# Pboot <- if (is.null(object$pboot)) { 
			# NULL
	   		# } else if (object$pboot < 1e-16) { 
			# "<1e-16"
	   		# } else if (object$pboot < 1e-4) { 
			# paste0("<1e", as.integer(ceiling(log(object$pboot, 10))))  
	   		# } else object$pboot
		
		# PtGKF <- if (is.null(object$ptGKF)) { 
		  # NULL
    		# } else if (object$ptGKF < 1e-16) { 
    		  # "<1e-16"
    		# } else if (object$ptGKF < 1e-4) { 
    		  # paste0("<1e", as.integer(ceiling(log(object$ptGKF, 10))))  
    		# } else object$ptGKF
    		
    	# if (object$scbtype == "normal") {
			# statresult <- data.frame(object$teststat, Pnorm)
			# names(statresult) <- c("supnorm", "p")
		# } else if (object$scbtype == "bootstrap") {
			# statresult <- data.frame(object$teststat, Pboot)
			# names(statresult) <- c("supnorm", "p")
		# } else if (object$scbtype == "tGKF") {
			# statresult <- data.frame(object$teststat, Pboot)
			# names(statresult) <- c("supnorm", "p")
		# } else if (object$scbtype == "all") { 
			# statresult <- data.frame(object$teststat, Pnorm, Pboot, PtGKF)
			# names(statresult) <- c("supnorm", "normal p", "bootstrap p", "tGKF p")	
		# }
		
		# if (length(object$model)) {		
		if (callfun == "scb.model")
			cat("Goodness-of-fit test\n") 	
			if (object$model == 0) {
				cat ("Model: zero mean function\n") 
			} else if (object$model == 1) {
				cat("Model: linear mean function\n")
			} else if (length(object$model) == 1) {
				cat("Model: polynomial mean function of degree <=", object$model,"\n")
			} else cat("Model: mean function in function space of dimension", ncol(object$model),"\n")	
			cat("Bandwidth:", round(object$bandwidth, 4), "\n")
			cat("Grid size:", object$gridsize, "\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", all = "normal and bootstrap"),"\n")
			if (object$scbtype != "no") {
				cat("Significance level:", object$level, "\n")
				cat("Replicates:", switch(object$scbtype, 
					normal = paste(object$nrep,"(normal)\n"), 
					bootstrap = paste(object$nboot,"(bootstrap)\n"), 
					all = paste(object$nrep,"(normal)", 
					object$nboot,"(bootstrap)\n")))
				if (scb.type != "no") {
					cat("Test statistic and p value(s)\n")
					print(statresult, right = FALSE, print.gap = 2L, digits = 4L,
					 	row.names = FALSE)
				}
			} 
	
		} 
		else {
			cat("Equality test for mean functions\n")
			cat("Bandwidths:", round(object$bandwidth,4),"\n")
			cat("Grid size:", object$gridsize, "\n")
			cat("SCB type:", switch(object$scbtype, no = "no SCB", normal = "normal", 
				bootstrap = "boostrap", tGKF = "tGKF", 
				all = "normal, bootstrap, and tGKF"),"\n")
			# if (object$scbtype != "no") {
			cat("Significance level:", object$level, "\n")
			cat("Replicates:", switch(object$scbtype, 
				normal = paste(object$nrep,"(normal)\n"), 
				bootstrap = paste(object$nboot,"(bootstrap)\n"), 
				all = paste(object$nrep,"(normal)", object$nboot,"(bootstrap)\n")))
				# if (scb.type != "no") {
			cat("Test statistic and p value(s)\n")
			print(statresult, right = FALSE, print.gap = 2L, digits = 4L,
						row.names = FALSE)
				# }
			# } 
		}
	}

}
