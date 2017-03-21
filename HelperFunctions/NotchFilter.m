function yFilt = NotchFilter( y, f, bandwidth, varargin )
  parser = inputParser();
  parser.addOptional( 'dT', 1 )
  okType = @(t) ismember( t, {'pass', 'reject'} );
  parser.addParameter( 'type', 'reject', okType )
  parser.addParameter( 'bothDirections', true )
  parser.parse( varargin{:} )
  options = parser.Results;
  
  dT = options.dT;
  f = f * dT;
  filterCoefs = getFilterCoefs( f, bandwidth, options.type );
  
  needTranspose = ~isrow( y );
  if needTranspose
    y = y';
  end
  numY = numel( y );
  p = polyfit( 1:numY, y, 1 );
  yTrend = polyval( p, 1:numY );
  y = y - yTrend;

  yFilt = symmetricFilter( y, filterCoefs, f );
  if options.bothDirections
    yFilt = flip( symmetricFilter( flip( yFilt ), filterCoefs, f ) );
  end
  
  if strcmpi( options.type, 'reject' )
    yFilt = yFilt + yTrend;
  end

  if needTranspose
    yFilt = yFilt';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function filterCoefs = getFilterCoefs( f, bandwidth, type )
  R = 1 - 3 * bandwidth;
  twoCosF = 2 * cos( 2 * pi * f );
  K = (1 - R * twoCosF + R * R) / (2 - twoCosF);

  switch lower( type )
    case 'reject'
      a0 = K;
      a1 = -K * twoCosF;
      a2 = K;
      b1 = R * twoCosF;
      b2 = -R*R;
    case 'pass'
      a0 = 1 - K;
      a1 = (K-R) * twoCosF;
      a2 = R * R - K;
      b1 = R * twoCosF;
      b2 = -R*R;
    otherwise
      error( 'Invalid band type: %s. Valid options are ''pass'' or ''reject''' )
  end
  filterCoefs = [a0, a1, a2, b1, b2];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function yFilt = symmetricFilter( y, filterCoefs, f )
  numY = numel( y );
  numPad = min( numY - 1, ceil( 2 / f ) );
  ySym = [y((numPad+1):-1:2), y];
  
  a0 = filterCoefs(1); a1 = filterCoefs(2); a2 = filterCoefs(3);
                       b1 = filterCoefs(4); b2 = filterCoefs(5);
  yFilt = zeros( 1, numY + numPad );
  yFilt(1) = a0 * ySym(1);
  yFilt(2) = a0 * ySym(2) ...
           + a1 * ySym(1) ...
           + b1 * yFilt(1);
  for n = 3:(numY+numPad)
    yFilt(n) = a0 * ySym(n) ...
             + a1 * ySym(n-1) ...
             + a2 * ySym(n-2) ...
             + b1 * yFilt(n-1) ...
             + b2 * yFilt(n-2);
  end
  yFilt(1:numPad) = [];  
end
