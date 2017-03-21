% [ySpectrum, f, yCorr] = Spectrum(y, dt, method)
function varargout = Spectrum(y, dT, varargin)
  parser = inputParser();
  parser.addParameter( 'removeTrend', true )
  parser.addParameter( 'plot', false )
  parser.addParameter( 'title', 'Spectrum' )
  parser.addParameter( 'xlabel', 'Frequency' )
  parser.addParameter( 'ylabel', 'Power' )
  parser.parse( varargin{:} )
  options = parser.Results;
  
  needTranspose = ~isrow( y );
  if needTranspose
    y = y';
  end
  if options.removeTrend
    % remove linear trend from y, unless requested not to
    y = removeTrend( y );
  end
  if isempty( dT )
    dT = 1;
  end
  
  % compute spectrum
  yCorr = autocorr( y );
  
  numY = numel( y );
  halfInd = ceil( numY / 2 );
  ySpectrum = fft( yCorr(halfInd:(halfInd + numY - 1)) );
  
  if options.plot || nargout > 1
    % need to compute the frequencies
    numSpec = numel( ySpectrum );
    halfInd = ceil( numSpec/2 );
    f = (0:(numSpec-1)) ./ (dT * numSpec);
    f(halfInd:end) = f(halfInd:end) - 1/(dT * numSpec);
    if options.plot
      yPower = sqrt( abs( ySpectrum ) );
      titleStr = options.title;
      fig = NamedFigure( titleStr, 'WindowStyle', 'docked' ); clf( fig )
      ax = axes( 'Parent', fig, 'OuterPosition', [0 0 1 1] );
      plot( ax, f(1:(halfInd-1)), yPower(1:(halfInd-1)) )
      title( ax, titleStr )
      ylabel( ax, options.ylabel )
      xlabel( ax, options.xlabel )
    end
  end
  
  if nargout == 0
    varargout = {};
  else
    varargout = {ySpectrum, f, yCorr};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove linear trend from y
function y = removeTrend( y )
  x = 1:numel( y );
  p = polyfit( x, y, 1 );
  y = y - polyval( p, x );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yCorr = autocorr( y )
  numY = numel( y ); numC = 2 * numY - 1;
  yFft = fft( y, numC );
  yCorr = ifft( yFft .* conj( yFft ), 'symmetric' );
  scale = numY - abs((1-numY):(numY-1));
  yCorr = yCorr ./ scale;
end