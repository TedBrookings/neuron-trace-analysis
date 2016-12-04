% [peak, sigmaMinus, sigmaPlus] = findPeak( x, sigmaCheckQuantile, outlierQuantile )
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
  if maxInd == numDensityPts
    halfInd = floor( numDensityPts / 2 );
    [maxDense, maxInd] = max( density(1:halfInd) );
  end
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
    sigmaPlus = (checkVal - peak) / numSigmaCheck;
    if maxInd > 1
      checkInd = [];
    end
    while isempty( checkInd )
      numSigmaCheck = sqrt( 2.0 ) * erfinv( sigmaCheckQuantile );
      sigmaDense = maxDense * exp( -numSigmaCheck^2 / sqrt(2) );
      checkInd = find( density(1:maxInd-1) <= sigmaDense, 1, 'last' );
      if isempty( checkInd )
        sigmaCheckQuantile = 0.9 * sigmaCheckQuantile;
      end
    end
    checkVal = vals(checkInd);
    sigmaMinus = abs(peak - checkVal) / numSigmaCheck;
  else
    sigmaPlus = []; sigmaMinus = [];
  end
  varargout = { peak, sigmaMinus, sigmaPlus };
end

