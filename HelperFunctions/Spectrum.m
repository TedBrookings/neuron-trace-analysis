% [powerSpectrum, f, yCorr] = Spectrum(y, dt, method)
function varargout = Spectrum(y, dT, varargin)
  parser = inputParser();
  parser.addParameter( 'removeTrend', true )
  parser.addParameter( 'plot', false )
  parser.addParameter( 'title', 'Spectrum' )
  parser.addParameter( 'xlabel', 'Frequency' )
  parser.addParameter( 'ylabel', 'Power' )
  parser.addParameter( 'yScale', 'log' )
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
  powerSpectrum = fft( yCorr(halfInd:(halfInd + numY - 1)) );
  halfInd = ceil( numY / 2 );
  powerSpectrum = sqrt( abs( powerSpectrum(1:halfInd) ) );
  
  if options.plot || nargout > 1
    % need to compute the frequencies
    f = (0:halfInd-1) ./ (dT * numY);
    if options.plot
      titleStr = options.title;
      fig = NamedFigure( titleStr, 'WindowStyle', 'docked' ); clf( fig )
      ax = axes( 'Parent', fig, 'OuterPosition', [0 0 1 1] );
      plot( ax, f, powerSpectrum )
      title( ax, titleStr )
      ylabel( ax, options.ylabel )
      xlabel( ax, options.xlabel )
      ax.YScale = options.yScale;
    end
  end
  
  if nargout == 0
    varargout = {};
  else
    varargout = {powerSpectrum, f, yCorr};
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
  yFft = fft( y ./ norm( y ), numC );
  yCorr = ifft( yFft .* conj( yFft ), 'symmetric' );
  scale = numY - abs((1-numY):(numY-1));
  yCorr = yCorr ./ scale;
end