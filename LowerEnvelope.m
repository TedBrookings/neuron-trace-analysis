% env = LowerEnvelope( signal, filtLength )
% Return lower envelope of signal
function env = LowerEnvelope( signal, filtLength, varargin )
  parser = inputParser();
  parser.addParameter('interpStyle', 'linear')
  parser.addParameter('plot', false)
  parser.addParameter('title', 'LowerEnvelope')
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  transpose = iscolumn( signal );
  if transpose
    signal = signal';
  end
  if ~isrow( signal )
    error( 'signal must be a 1D array' )
  end
  
  nanSignal = isnan( signal );
  if any( nanSignal )
    % break signal up into segments of non-nan data, and get lower envelope
    % of good segments
    env = nan( size(signal) );
    starts = find( ~nanSignal & [ true, nanSignal(1:end-1) ] );
    stops = find( ~nanSignal & [ nanSignal(2:end), true ] );
    for n = 1:numel( starts )
      start = start(n); stop = stops(n);
      env(start:stop) = getLowerEnvelope( signal(start:stop), ...
                                          filtLength, options );
    end
  else
    env = getLowerEnvelope( signal, filtLength, options );
  end
  
  if options.plot
    plotLowerEnvelope( signal, env, options )
  end
  
  if transpose
    env = env';
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the lower envelop of signal with no NaN data
function env = getLowerEnvelope( signal, filtLength, options )
  env = imerode( signal, ones( 1, filtLength ) );
  
  if strcmpi( options.interpStyle, 'none' )
    return % no interpolation
  end

  % find points where the envelope is >= signal
  touchPoints = getTouchPoints( signal, env );
  % interpolate the envelope between touchPoints
  env = interp1( touchPoints, signal(touchPoints), 1:numel( signal ), ...
                 options.interpStyle );
  % find height of env above signal (indicates need for refinement)
  err = env - signal;
  while any( err > 0 )
    [maxErr, errInds] = arrayfun( @(i1,i2) checkMax( err(i1:i2) ), ...
                                touchPoints(1:end-1), touchPoints(2:end) );
    errInds = errInds + touchPoints(1:end-1) - 1;
    errInds( maxErr <= 0 ) = [];
    touchPoints = sort( [touchPoints, errInds] );
    env = interp1( touchPoints, signal(touchPoints), 1:numel( signal ), ...
                   options.interpStyle );
    % find height of env above signal (indicates need for refinement)
    err = env - signal;
  end
  
  function [maxVal, maxInd] = checkMax( vals )
    if isempty(vals)
      maxVal = zeros( 'like', vals ); maxInd = 0;
    else
      [maxVal, maxInd] = max( vals );
    end
  end
  %{
  converged = false;
  while ~converged
    % interpolate the envelope between touchPoints
    env = interp1( touchPoints, signal(touchPoints), 1:numel( signal ), ...
                   options.interpStyle );
    % save old, get new touchPoints
    oldTouchPoints = touchPoints;
    touchPoints = getTouchPoints( signal, env );
    % loop until the touchPoints don't change
    converged = isequal( touchPoints, oldTouchPoints );
  end
%}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find points where the envelope is >= signal.
function touchPoints = getTouchPoints( signal, env )
  touchPoints = find( env >= signal );
  % insist that the end-points be included
  if touchPoints(1) > 1
    touchPoints = [ 1 touchPoints ];
  end
  if touchPoints(end) < numel( signal )
    touchPoints = [ touchPoints numel( signal ) ];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% visualize the results of the process
function plotLowerEnvelope( signal, env, options )
  fig = NamedFigure( options.title ); fig.WindowStyle = 'docked';
  clf( fig )
  
  ax1 = subplot( 2, 1, 1, 'Parent', fig );
  plot( ax1, signal, 'b-' )
  hold( ax1, 'on' )
  plot( ax1, env, 'r-', 'LineWidth', 2 )
  ax1.XTick = [];
  axis( ax1, 'tight' )
  lg1 = legend( ax1, { 'signal', 'envelope' }, 'Location', 'NorthOutside', ...
          'Orientation', 'Horizontal', 'Units', 'normalized', ...
          'FontSize', 8 );
  setPos( ax1, lg1, [ 0.05, 0.53, 0.88, 0.5/1.2 ] )
  
  ax2 = subplot( 2, 1, 2, 'Parent', fig );
  plot( ax2, signal - env, 'b-' )
  axis( ax2, 'tight' )
  lg2 = legend( ax2, { 'signal - envelope' }, 'Location', 'NorthOutside', ...
               'Orientation', 'Horizontal', 'Units', 'normalized', ...
               'FontSize', 8 );
  setPos( ax2, lg2, [ 0.05, 0.07, 0.88, 0.5/1.2 ] )
  
  linkaxes( [ax1, ax2], 'x' )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function setPos( ax, lg, pos )
  ax.Units = 'normalized';
  ax.FontUnits = 'normalized';
  ax.FontSize = 0.05;
  ax.Position = pos;
  
  lg.Position(2) = sum( ax.Position([2 4]) );
end