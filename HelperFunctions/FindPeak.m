% [peak, sigma] = findPeak( x, sigmaCheckQuantile, outlierQuantile )
% find the peak of list of data points x
% NOTE - this function assumes x is a sorted row of finite values
%   INPUTS:
%      -x: data values
%      -sigmaCheckQuantile: what quantile to look at to determine sigma
%      -outlierQuantile: ignore data this quantile difference from 0 or 1
%                        (defaults to 0.05)
function varargout = FindPeak( x, sigmaCheckQuantile, outlierQuantile )
  if nargin < 3
    outlierQuantile = 0.05;
  end
  numPts = numel( x );
  i1 = 1 + round( (numPts - 1) * outlierQuantile );
  i2 = 1 + round( (numPts - 1) * (1 - outlierQuantile) );
  numDensityPts = max( 100, round( sqrt( numPts ) ) );
  vals = linspace( x(i1), x(i2), numDensityPts );
  density = ksdensity( x, vals );
  [maxDense, maxInd] = max( density );
  peak = vals(maxInd);

  if nargout > 1
    checkInd = [];
    while isempty( checkInd )
      numSigmaCheck = sqrt( 2.0 ) * erfinv( sigmaCheckQuantile );
      sigmaDense = maxDense * exp( -numSigmaCheck^2 / sqrt(2) );
      checkInd = maxInd + find( density(maxInd+1:end) <= sigmaDense, 1 );
      if isempty( checkInd )
        sigmaCheckQuantile = 0.9 * sigmaCheckQuantile;
      end
    end
    checkVal = vals(checkInd);
    sigma = peak + (checkVal - peak) / numSigmaCheck;
  else
    sigma = [];
  end
  varargout = { peak, sigma };
end

