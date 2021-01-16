# Returns the expected Euler characteristic curve for a t-field if its
# Lipschitz-Killing-curvatures and degrees of freedom are provided.
#
# @param LKC Vector estimated LKCs of the UGRFs of the field.
# @param df Numeric degrees of freedom of the field.
# @return a function computing the EEC curve
# @export
EEC <- function( LKC, df ){
  # Approximate the tail distribution using the GKF and EECH
  EECf <- function( x ){
              1   * ( 1 - stats::pt( x, df = df ) ) +
              LKC * ( 2*pi )^(-1) * ( 1 + x^2 / df )^( -( df - 1 ) / 2 )
  }
  # Return the vectorized function
  Vectorize( EECf )
}

EEC_threshold <- function( LKC,
                           alpha    = 0.05,
                           df,
                           interval = c( 0, 100 )
                           ){
  # Get the EEC curve and subtract alpha
  EEC_function <- EEC( LKC, df )
  tailProb     <- function( u ){ EEC_function( u ) - alpha }

  # return the approximation of the quantile from the tGKF
  list( q = stats::uniroot( tailProb, interval )$root, EEC = EEC_function )
}

# Utility function: spline interpolation + numerical derivative
splineDeriv <- function(y, x, h = NULL)
{
	if (is.null(h)) h <- .Machine$double.eps^(1/3)
	yout <- stats::spline(x, y, method="natural", xout=c(x-h,x+h))$y
	dim(yout) <- c(length(x),2)
	return((yout[,2]-yout[,1])/(2*h))
}


# Estimates the LKCs from normed residuals
LKCest <- function( R, x = NULL, Q = NULL, xQ = NULL ){
  # Get default coordinate values
  if( is.null(x) ){
    # x <- seq( 0, 1, length.out = dimR[1] )
    x <- seq( 0, 1, length.out = dim(R)[1] )
  }
  
    
  #---------------------------------------------------------------------------
  # Compute L1
  #---------------------------------------------------------------------------
  # Observed location on the connected component
  # Get the derivative of the residuals
  # dR <- apply( R, 2,
               # FUN = function( yy ){
                     # # Interpolate the function
                     # fn <- stats::splinefun( x = x,
                                             # y = yy,
                                             # method = "natural" )
                     # # Obtain the derivative of the function
                     # pracma::fderiv( f = fn,
                                     # x = x,
                                     # n = 1,
                                     # method = "central")
                   # } )                  
  dR <- apply(R, 2, splineDeriv, x = x)

  # get standard deviation of derivative
  dR.sd <- apply( dR , 1, stats::sd )
  
  if( is.null(Q) ){
    # integrate the standard deviation of the derivative using trapezoid rule to
    # get the LKC
    L <- sum( diff(x) / 2 * ( dR.sd[ 1:( length(dR.sd) - 1 ) ] + dR.sd[ 2:length(dR.sd) ] ) )
    
    return( L )
    
  }else{
    #---------------------------------------------------------------------------
    # Compute L1
    #---------------------------------------------------------------------------
    # observed location on the connected component
    # get the derivative of the residuals
    # dQ <- apply( Q, 2,
                 # FUN = function( yy ){
                   # # Interpolate the function
                   # fn <- stats::splinefun( x = xQ,
                                           # y = yy,
                                           # method = "natural" )
                   # # Obtain the derivative of the function
                   # pracma::fderiv( f = fn,
                                   # x = x,
                                   # n = 1,
                                   # method = "central")
                 # } )
    dQ <- apply(Q, 2, splineDeriv, x = xQ)
             
    # get standard deviation of derivative
    dQ.sd <- apply( dQ , 1, stats::sd )
    
    # integrate the standard deviation of the derivative using trapezoid rule to
    # get the LKC
    dvar = sqrt( dR.sd^2 + dQ.sd^2 )
    L <- sum( diff(x) / 2 * ( dvar[ 1:( length(dR.sd) - 1 ) ] +
                              dvar[ 2:length(dR.sd) ] ) )
    
    return( L )
    
    
  }
}
