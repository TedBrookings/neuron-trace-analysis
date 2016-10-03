function filterFunc = GetFilterFunction( filterScales, varargin )
  parser = inputParser();
  parser.addParameter( 'filterDim', [] )
  parser.addParameter( 'filterGaussianSigmaRatio', 6.0 )
  parser.addParameter( 'filterNumSigma', 2.0 )
  parser.addParameter( 'filterVector', [] )

  parser.parse( varargin{:} )
  options = parser.Results;
  
  if isempty( options.filterVector )
    filterVector = getGaussianFilter( filterScales(1), options );
    for n = 2:numel( filterScales )
      f_n = getGaussianFilter( filterScales(n), options );
      filterVector = addFilters( filterVector, f_n );
    end
  else
    filterVector = options.filterVector;
  end
  
  if ~iscolumn( filterVector )
    filterVector = filterVector';
  end

  filterFunc = @(y) applyFilter( y, filterVector, options.filterDim );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function gaussianFilter = getGaussianFilter( filterScale, options )
  filterLength = options.filterGaussianSigmaRatio * abs( filterScale );
  halfLength = round( (filterLength - 1) / 2 );
  if halfLength < 1
    gaussianFilter = 1;
  else
    filterLength = 1 + 2 * halfLength;
    x = linspace( -options.filterNumSigma, options.filterNumSigma, ...
                  filterLength );
    gaussianFilter = exp( -x.^2 ./ 2 );
    gaussianFilter = gaussianFilter ./ sum( gaussianFilter );
  end
  
  if filterScale < 0
    gaussianFilter = -gaussianFilter;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filterVec = addFilters( filterVec1, filterVec2 )
  num1 = numel( filterVec1 ); num2 = numel( filterVec2 );
  if num1 > num2
    filterVec = filterVec1;
    deltaInd = (num1 - num2) / 2;
    i1 = 1 + deltaInd;
    i2 = num1 - deltaInd;
    filterVec(i1:i2) = filterVec(i1:i2) + filterVec2;
  else
    filterVec = filterVec2;
    deltaInd = (num2 - num1) / 2;
    i1 = 1 + deltaInd;
    i2 = num2 - deltaInd;
    filterVec(i1:i2) = filterVec(i1:i2) + filterVec1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% given a signal y, apply convolution filter represented by filterVec along
% dimension filterDim
function y = applyFilter( y, filterVec, filterDim )
  if isempty( filterDim )
    filterDim = find( size( y ) > 1, 1 );
    if isempty( filterDim )
      return
    end
  end
  switch filterDim
    case 1, needTranspose = false;
    case 2, needTranspose = true;
    otherwise, error( 'Filter only works for vectors or 2D matrices' )
  end
  if needTranspose, y = y'; end
  needConvert = ~isa( y, 'double' );
  if needConvert
    yClass = class( y );
    y = double( y );
  end
  if ~iscolumn( filterVec ), filterVec = filterVec'; end
  
  % use convolution to filter signal, but pad first so that only valid part
  % of convolution is used
  y = conv2( padSignal( y, filterVec ), filterVec, 'valid' );

  if needTranspose, y = y'; end
  if needConvert, y = feval( yClass, y ); end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% given a signal y and a convolution filter represented by filterVec, pad y
% symmetrically so that the desired filtered signal is the "valid" part of
% the convolution
function y = padSignal( y, filterVec )
  filtLen = numel( filterVec );
  halfLen = (filtLen - 1) / 2;
  yLen = size( y, 1 ); originalYLen = yLen;
  finalYLen = originalYLen + 2 * halfLen;
  while yLen < finalYLen
    if 3 * yLen - 2 <= finalYLen
      % pad y fully with itself. Will complete pad if above is an equality
      pad = flip( y, 1 );
      y = [ pad(1:end-1) ; y ; pad(2:end) ]; %#ok<AGROW>
    else
      % pad y partially: this should complete pad
      padLen = (finalYLen - yLen) / 2;
      frontPadLen = floor( padLen );
      backPadLen = ceil( padLen );
      frontPad = y(1+frontPadLen:-1:2,:);
      backPad = y(end-1:-1:end-backPadLen,:);
      y = [ frontPad ; y ; backPad ]; %#ok<AGROW>
    end
    yLen = size( y, 1 );
  end
end