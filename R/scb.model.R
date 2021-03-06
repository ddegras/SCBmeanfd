scb.model <- function(x, y, model, bandwidth, level = .05, degree = 1,
	scbtype = c("normal","bootstrap","tGKF","all"), gridsize = length(x), 
	keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no", "multicore", "snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	stopifnot(all(!is.na(x) & !is.na(y)))
	stopifnot(length(x) == ncol(y)) 
	stopifnot(NROW(model) %in% c(1,length(x))) 

	n <- nrow(y)
	if (length(model) == 1 && model == 0) {
		par.res <- t(y)
		par.mu.hat <- rep(0, gridsize)
	} 
	else if (length(model) == 1 && model > 0) {
		par.y.hat <- lm( t(y) ~ poly(x, degree = model))
		par.res <- residuals(par.y.hat)
		par.mu.hat <- rowMeans(fitted(par.y.hat))
	} 
	else {
		par.y.hat <- lm( t(y) ~ model - 1) 
		par.res <- residuals(par.y.hat)		    
		par.mu.hat <- rowMeans(fitted(par.y.hat))	
	}
	smooth.par.mu.hat <- locpoly(x, par.mu.hat, degree = degree, 
		bandwidth = bandwidth, gridsize = gridsize)$y
	  nonpar.mu.hat <- locpoly(x, colMeans(y), degree = degree, 
		bandwidth = bandwidth, gridsize = gridsize)$y
	r <- apply(par.res, 2, function(z) locpoly(x, z, degree = degree, 
		bandwidth = bandwidth, gridsize = gridsize)$y)

	rbar <- rowMeans(r)
	r <- r - rbar
	sigma.hat <- sqrt(rowSums(r^2) / (n-1))
	se <- sigma.hat / sqrt(n)
	test.stat <- max(abs(rbar/se))
	p.norm = p.boot = p.tGKF = NULL
	q.norm = q.boot = q.tGKF = NULL
	lb.norm = ub.norm = lb.boot = ub.boot = lb.tGKF = ub.tGKF = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","all")) {
		svd.r <- svd(r / (sigma.hat * sqrt(n-1)), nv = 0)
		ncomp <- which.max(cumsum(svd.r$d^2) > .99 * sum(svd.r$d^2))
		vars <- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- t(t(svd.r$u[,1:ncomp]) * svd.r$d[1:ncomp])
		supnorm <- apply(abs(M %*% vars), 2, max)
		p.norm <- 1 - ecdf(supnorm)(test.stat)
		q.norm <- as.numeric(quantile(supnorm,1-level)) 
		lb.norm <- nonpar.mu.hat - q.norm * se
		ub.norm <- nonpar.mu.hat + q.norm * se	
	}

	if (scbtype %in% c("bootstrap","all")) {
		boot.stat <- function(mat,ix) {
			e.boot <- colMeans(mat[ix,])
			sqrt((n-1) * max((e.boot^2)/(colSums(mat[ix,]^2)-n*e.boot^2)) * n)
		}
		supnorm <- boot(t(r), boot.stat, nboot, 
			parallel = parallel, ncpus = ncpus, cl = cl)$t 
		p.boot <- 1 - ecdf(supnorm)(test.stat)
		q.boot <- as.numeric(quantile(supnorm,1-level))
		lb.boot <- nonpar.mu.hat - q.boot * se
		ub.boot <- nonpar.mu.hat + q.boot * se	
	}
	
	if (scbtype %in% c("tGKF","all")) {
	  # Estimate the LKCs
	  xgrid <- seq(min(x), max(x), len=gridsize)
	  L = LKCest( R = r / sigma.hat, x = xgrid )
	  # Get the tGKF threshold
	  q.tGKF <- EEC_threshold( L,
	                           alpha    = ( 1 - level ) * 0.5,
	                           df = n - 1,
	                           interval = c( 0, 100 )
	  )
	  p.tGKF <- min(q.tGKF$EEC(test.stat),1)
	  q.tGKF <- q.tGKF$q
	  
	  # Get the simultaneous confidence bands
	  lb.tGKF <- nonpar.mu.hat - q.tGKF * se
	  ub.tGKF <- nonpar.mu.hat + q.tGKF * se	
	}

	result <- list( x = x, y = if(keep.y) y else NULL, call = caLL, model = model, 
		par = smooth.par.mu.hat, nonpar = nonpar.mu.hat, bandwidth = bandwidth, 
		degree = degree, level = level, scbtype = scbtype, teststat = test.stat,
		pnorm = p.norm, pboot = p.boot,  ptGKF = p.tGKF, 
		qnorm = q.norm, qboot = q.boot, qtGKF = q.tGKF, 
		normscb = cbind(lb.norm, ub.norm), bootscb = cbind(lb.boot, ub.boot), 
		tGKFscb = cbind(lb.tGKF, ub.tGKF), gridsize = gridsize, nrep = nrep, nboot = nboot )

	class(result) <- "SCBand"
	return(result)	
		
}
