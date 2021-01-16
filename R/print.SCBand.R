print.SCBand <- function(object, ...) 
{
	stopifnot(inherits(object,"SCBand"))
	callfun <- as.character(object$call[1])
	if (callfun == "scb.mean") {
		cat("\nMean function estimation\n")
		cat("Bandwidth:", round(object$bandwidth,4), "\n")
		cat("SCB type:", switch(object$scbtype, normal = "normal", 
			bootstrap = "boostrap", tGKF = "tGKF", 
			all = "normal, bootstrap, tGKF"),"\n")
		cat("Confidence level:", object$level, "\n")
		cat("Quantile used for SCB:\n")		
		statresult <- switch(object$scbtype,
			normal = data.frame(normal = object$qnorm),
			boot = data.frame(bootstrap = object$qboot),
			tGKF = data.frame(tGKF = object$qtGKF),
			all = data.frame(normal = object$qnorm, 
				bootstrap = object$qboot, tGKF = object$qtGKF))
		print(statresult, print.gap = 2L, right = FALSE, digits = 4L, 
			row.names = FALSE)		
	} else {		
		pfun <- function(p) {
			if (is.null(p)) return(NA)
			if (p < 1e-16) return("< 1e-16")
			if (p < 1e-4) return(paste0("<1e", 
				as.integer(ceiling(log(p, 10)))))
			return(p) 
		}
		statresult <- switch(object$scbtype,
			normal = data.frame(object$teststat, pfun(object$pnorm)),
			bootstrap = data.frame(object$teststat, pfun(object$pboot)),
			tGKF = data.frame(object$teststat, pfun(object$ptGKF)),
			all = data.frame(object$teststat, pfun(object$pnorm), 
				pfun(object$pboot), pfun(object$ptGKF)))
		names(statresult) <- switch(object$scbtype,
			normal = c("stat", "normal-p"),
			bootstrap = c("stat", "bootstrap-p"),
			tGKF = c("stat", "tGKF-p"),
			all = c("stat", "normal-p", "bootstrap-p", "tGKF-p"))
		
		if (callfun == "scb.model") {
			cat("\nGoodness-of-fit test\n") 	
			cat ("Model for the mean function: ")
			if (length(object$model) == 1) {
				if (object$model == 0) { cat ("zero\n") 
				} else if (object$model == 1) { cat("linear\n")
				} else cat("polynomial of degree <=", object$model,"\n")
			} else cat("function space of dimension", ncol(object$model),"\n")
			cat("Bandwidth:", round(object$bandwidth, 4), "\n")
			cat("SCB type:", switch(object$scbtype, normal = "normal", 
				bootstrap = "boostrap", tGKF = "tGKF", 
				all = "normal, bootstrap, tGKF"),"\n")
			cat("Significance level:", object$level, "\n")
			cat("Test statistic and p-value\n")
			print(statresult, right = FALSE, print.gap = 2L, digits = 4L, 
				row.names = FALSE)
		} else {
			cat("\nEquality test for mean functions\n\n")
			cat("Bandwidths:", round(object$bandwidth,4),"\n")
			cat("SCB type:", switch(object$scbtype, normal = "normal", 
				bootstrap = "boostrap", tGKF = "tGKF", 
				all = "normal, bootstrap, tGKF"),"\n")
			cat("Significance level:", object$level, "\n")
			cat("Test statistic and p-value\n")
			print(statresult, right = FALSE, print.gap = 2L, digits = 4L, 
				row.names = FALSE)
		} 
	}
}
