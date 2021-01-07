cv.select <- function (x, y, degree = 1, interval = NULL, 
	gridsize = length(x), ...) 
{
    if (is.null(interval)) {
        rangex <- diff(range(x))
        meshx <- max(diff(sort(x)))
        interval <- c(ifelse(degree < 2, meshx/2, meshx), rangex/2)
    }
    
    ## Pilot grid search
    h.gridsize <- 11
    h.grid <- seq(interval[1], interval[2], len=h.gridsize)
    score.grid <- sapply(h.grid, cv.score, x=x, y=y, 
    		degree = degree, gridsize = gridsize, ...)
    idx <- which.min(score.grid)
    
    ## Refined line search
    lo <- h.grid[max(idx-1,1)]
    hi <- h.grid[min(idx+1,h.gridsize)]
    optimize(cv.score, c(lo,hi), x = x, y = y, 
		degree = degree, gridsize = gridsize, ...)$minimum
}
