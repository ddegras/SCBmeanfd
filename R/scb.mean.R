scb.mean <- function(x, y, bandwidth, level = .95, degree = 1, 
	scbtype = c("normal","bootstrap","both","tGKF","no"), gridsize = length(x), 
	keep.y = TRUE, nrep = 2e4, nboot = 5e3, parallel = c("no","multicore","snow"), 
	ncpus = getOption("boot.ncpus",1L), cl = NULL)
{
	caLL <- match.call()
	stopifnot(is.matrix(y))
	stopifnot(all(!is.na(x) & !is.na(y)))
	stopifnot(length(x) == ncol(y)) 
	n <- nrow(y)
	N <- ncol(y)	
	y.hat <- apply(y, 1, function(z) locpoly(x, z, degree = degree, 
		bandwidth = bandwidth, gridsize = gridsize)$y)
	mu.hat <- rowMeans(y.hat)
	sigma.hat <- apply(y.hat, 1, sd)
	se <- sigma.hat / sqrt(n)
	r <- y.hat - mu.hat  	
	lb.norm = ub.norm = lb.boot = ub.boot = q.norm = q.boot = NULL
	scbtype <- match.arg(scbtype)	

	if (scbtype %in% c("normal","both")) {
		svd.r <- svd(r / sigma.hat / sqrt(n-1), nv = 0)
		ncomp <- which.max(cumsum(svd.r$d^2) > .99 * sum(svd.r$d^2))
		vars <- matrix(rnorm(ncomp * nrep), ncomp, nrep)
		M <- svd.r$u[,1:ncomp] %*% diag(svd.r$d[1:ncomp],ncomp,ncomp)
		supnorm <- apply(abs(M %*% vars), 2, max)
		q.norm <- as.numeric(quantile(supnorm,level)) 
		lb.norm <- mu.hat - q.norm * se
		ub.norm <- mu.hat + q.norm * se	
	}
	if (scbtype %in% c("bootstrap","both")) {
		boot.stat <- function(mat,ix) {
			mu.boot <- colMeans(mat[ix,])
			sqrt(max(((n-1)*mu.boot^2)/(colSums(mat[ix,]^2)-n*mu.boot^2)) * n)
		}
		supnorm   <- boot(t(r), boot.stat, nboot, parallel = parallel, 
			ncpus = ncpus, cl = cl)$t
		q.boot 	  <- as.numeric(quantile(supnorm,level)) 
		lb.boot   <- mu.hat - q.boot * se
		ub.boot   <- mu.hat + q.boot * se	
	}
	if (scbtype %in% c("tGKF","both")) {
	  # Estimate the LKCs
	  L = LKCest( R = r / se / sqrt(n), x = x )
	  # Get the tGKF threshold
	  q.tGKF <- EEC_threshold <- function( LKC,
	                                       alpha    = ( 1 - level ) * 0.5,
	                                       df = n - 1,
	                                       interval = c( 0, 100 )
	  ) 
	  lb.tGKF <- mu.hat - q.tGKF * se
	  ub.tGKF <- mu.hat + q.tGKF * se	
	}

	result <- list( x = x, y = if(keep.y) y else NULL, call = caLL, model = NULL,
					par      = NULL, nonpar = mu.hat, bandwidth = bandwidth,
					degree   = degree, level = level, scbtype = scbtype,
					teststat = NULL, pnorm = NULL, pboot = NULL, 
					qnorm    = q.norm, qboot = q.boot, qtGKF = q.tGKF,
					normscb  = cbind(lb.norm, ub.norm), 
					bootscb  = cbind(lb.boot, ub.boot),
					tGKFscb  = cbind(lb.tGKF, ub.tGKF),
					gridsize = gridsize, nrep = nrep, nboot = nboot )
	class(result) <- "SCBand"
	return(result)
}
