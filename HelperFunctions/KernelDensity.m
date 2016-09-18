%[density, samplePoints] = KernelDensity( data, samplePoints, varargin )
% Calculte kernel density (smooth estimate of density of samples) at samplePoints.
% INPUTS:
%   data: vector of input data. Non-finite data is ignored
%   samplePoints: points to estimate kernel density. If omitted, a sensible range
%                 is chosen based on the data
% OPTIONS
%    smoothingScale: (default=0.2) sets scale of gaussian smoothing (combined
%                    with number of data points and inter-sample-interval).
%                    Larger numbers will produce a smoother kernel estimate
% OUTPUTS
%   density: estimate of density at samplePoints
%   samplePoints: points where density is estimated
% Copyright 2016 Ted Brookings (ted.brookings@googlemailcom)
function varargout = KernelDensity( data, samplePoints, varargin )
  % process input options
  if nargin < 2
    samplePoints = [];
  end
  parser = inputParser();
  parser.addParameter( 'smoothingScale', 0.15 )
  parser.addParameter( 'adaptive', false )
  parser.addParameter( 'adaptiveNumSigma', 10.0 )
  
  parser.parse( varargin{:} )
  options = parser.Results;

  % process data to be sorted, finite, row-shaped
  data(~isfinite( data )) = [];
  if numel( data ) < 2
    density = zeros( size( samplePoints ) );
    varargout = { density, samplePoints };
    return
  end
  data = sort( data );
  dataTranspose = ~isrow( data );
  if dataTranspose, data = data'; end
  if isempty( samplePoints )
    pointsTranspose = dataTranspose;
  else
    pointsTranspose = ~isrow( samplePoints );
    if pointsTranspose, samplePoints = samplePoints'; end
  end
  
  kernelSigma = getKernelSigma( data, options );

  if isempty( samplePoints )
    samplePoints = getSamplePoints( data, kernelSigma );
  end

  % calculate density
  density = calculateDensity( data, samplePoints, kernelSigma, options );

  if dataTranspose, density = density'; end
  if pointsTranspose, samplePoints = samplePoints'; end
  varargout = { density, samplePoints };
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kernelSigma = getKernelSigma( data, options )
  % average Inter-Sample Interval
  averageIsi = median( diff( data ) );
  range = data(end) - data(1); numData = numel( data );
  numIsi = numData - 1;
  if ~( averageIsi * numIsi / range > 1.0e-6 )
    averageIsi = range / numIsi;
  end
  kernelSigma = options.smoothingScale * averageIsi * numData^0.8;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function samplePoints = getSamplePoints( data, kernelSigma )
  % average Inter-Sample Interval
  left = data(1); right = data(end);
  numSamplePoints = round( 3 * (right - left) / kernelSigma );
  numSamplePoints = max( min( numSamplePoints, 10000 ), 100 );
  samplePoints = linspace( left, right, numSamplePoints );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate density
function density = calculateDensity( data, samplePoints, kernelSigma, options )
  if options.adaptive
    nonadaptiveSigma = kernelSigma;
    neighborDist = kernelSigma * options.adaptiveNumSigma;
    numData = numel( data );
    for n = numData:-1:1
      [i1, i2] = bracketData( data, n, numData, neighborDist );
      kernelSigma(n) = getAdaptiveKernelSigma( data, i1, i2, nonadaptiveSigma, options );
      alpha(n) = -1.0 / (2 * kernelSigma(n)^2);
      gaussNorm(n) = sqrt( 2 * pi ) * numData * kernelSigma(n);
    end
    density = arrayfun( @(x) sum( exp( alpha .* (data - x).^2 ) ./ gaussNorm ), ...
                        samplePoints );
  else
    alpha = -1.0 / (2 * kernelSigma^2);
    gaussNorm = sqrt( 2 * pi ) * numel( data ) * kernelSigma;
    density = arrayfun( @(x) sum( exp( alpha .* (data - x).^2 ) ) / gaussNorm, ...
                        samplePoints );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find indices to data within neighborDist of data(n)
function [i1, i2] = bracketData( data, n, numData, neighborDist )
  dn = data(n); i1 = n; i2 = n;
  while i2 < numData && data(i2+1) - dn <= neighborDist
    i2 = i2 + 1;
  end
  while i1 > 1 && dn - data(i1-1) <= neighborDist
    i1 = i1 - 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function kernelSigma = getAdaptiveKernelSigma( data, i1, i2, nonadaptiveSigma, options )
  % average Inter-Sample Interval
  minAdaptiveNeighbors = 3;
  if i2 - i1 < minAdaptiveNeighbors
    kernelSigma = nonadaptiveSigma;
  else
    iData = data(i1:i2);
    averageIsi = median( diff( iData ) );
    if averageIsi * nonadaptiveSigma > 1.0e-6
      numData = numel( data );
      kernelSigma = options.smoothingScale * averageIsi * numData^0.8;
    else
      kernelSigma = nonadaptiveSigma;
    end
  end
end


