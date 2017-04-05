% fig = PlotSpikes( t, v, spike, burst, varargin )
% plot spike info (and optionally bursts) on top of voltage trace
%  INPUTS:
%   -t: 1xtraceLen time array (ms) or sampleInterval (ms)
%   -v: voltage trace array (mV)
%   -spike: structure with spike information, as calculated by GetSpikes.m
%   -burst: structure with burst information, as caculated by AnalyzeBurst.m
%           (if omitted or empty, no bursts are plotted)
function varargout = PlotSpikes( dt, v, spike, burst, varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;

  parser.addParameter( 'plotSubject', '' )
  parser.addParameter( 'timesOnly', false )
  parser.addParameter( 'spikeShapeMarkerSize', 6 )
  parser.addParameter( 'labelFontSize', 18 )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  options.plotBurst = exist( 'burst', 'var' ) && ~isempty( burst ) ...
            && ~isequal( burst, struct() );

  % get time vector (and convert units to seconds)
  [t, dt] = getTimeVec( dt, v );
  
  % create figure, axes, and title for plotting
  [fig, ax, titleStr] = createFig( options );
  
  % plot bursts if there are any detected
  if options.plotBurst
    legendEntries = plotBurst( ax, burst, v );
  else
    legendEntries = {};
  end
  
  % draw the voltage trace in white:
  voltageLine = plot( ax, t, v, 'w-' );
  % Exclude line from legend
  voltageLine.Annotation.LegendInformation.IconDisplayStyle = 'off';
  
  if options.timesOnly
    % plot spikes, but only show the spike times
    spikeLegendEntries = plotSpikesTimesOnly( ax, t, v, spike );
  else
    % plot spikes with spike shape info
    spikeLegendEntries = plotSpikesShapeInfo( ax, t, v, spike, options );
  end  
  legendEntries = [ legendEntries, spikeLegendEntries ];
  
  if ~isempty( legendEntries )
    legend( ax, legendEntries{:}, 'Location', 'SouthOutside', ...
            'Orientation', 'Horizontal' )
  end
  
  % add title and labels to axes
  addLabels( ax, titleStr, options )
  hold( ax, 'off' );
  if nargout > 0
    varargout = { fig, ax };
  else
    varargout = {};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get time vector (and convert units to seconds)
function [t, dt] = getTimeVec( dt, v )
  if numel( dt ) == 1
    % passed dt
    dt = dt / 1000; % convert to sec
    tFinal = dt * (numel( v ) - 1);
    t = 0:dt:tFinal;
  else
    % passed t (as dt)
    t = 0.001 * dt;
    if t(1) ~= 0
      t = t - t(1);
    end
    dt = (t(end) - t(1)) / (numel( v ) - 1);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create figure, axes, and title for plotting
function [fig, ax, titleStr] = createFig( options )
  if options.plotBurst
    % using burst information
    baseTitle = 'Spike/Burst Detection';
  else
    % only using spike information
    baseTitle = 'Spike Detection';
  end
  % make the title
  if ischar( options.plotSubject ) && ~isempty( options.plotSubject )
    titleStr = [options.plotSubject, ': ', baseTitle];
  else
    titleStr = baseTitle;
  end
  
  % create a new figure, name it, and make it ready for plotting
  fig = NamedFigure( titleStr ); clf( fig )
  fig.WindowStyle = 'docked';
  whitebg(fig, 'k');  % make the background black
  ax = subplot( 1,1,1, 'Parent', fig ); hold( ax, 'on' )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot bursts if there are any detected
function legendEntries = plotBurst( ax, burst, v )
  % get burst times
  burstTimes = burst.Times ./ 1000;
  burstLen = burst.Durations.list ./ 1000;
  % first draw blue rectangle to signify burst times
  numBursts = numel( burstTimes );
  if numBursts > 0
    % define upper and lower voltages of lines indicating burst times
    [top, bottom] = getAnnotationTopAndBottom( v );
    
    burstGroup = hggroup( 'parent', ax );
    for n = 1:numBursts
      tLow = burstTimes(n);
      tHigh = tLow + burstLen(n);
      fill( [tLow, tHigh, tHigh, tLow], [bottom, bottom, top, top], ...
            'b', 'Parent', burstGroup );
    end
    burstGroup.Annotation.LegendInformation.IconDisplayStyle = 'on';
    legendEntries = {'Burst'};
  else
    legendEntries = {};
  end
end  
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define upper and lower voltages of lines indicating burst times
function [top, bottom] = getAnnotationTopAndBottom( v )
  top = max( v );
  bottom = min( v );
  delta = 0.1 * (top - bottom);
  bottom = bottom - delta;
  top = top + delta;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spikes, but only show the spike times
function legendEntries = plotSpikesTimesOnly( ax, t, v, spike )
  % no shape information, just mark spikes
  if isempty( spike.ind )
    legendEntries = {};
  else
    plot( ax, t(spike.ind), v(spike.ind), 'ro', 'MarkerFaceColor', 'r' )
    legendEntries = {'spike peak'};
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot spikes with spike shape info
function legendEntries = plotSpikesShapeInfo( ax, t, v, spike, options )
  if isempty( spike.maxV.t )
    legendEntries = {};
    return
  end
  markerSize = options.spikeShapeMarkerSize;
  plot( ax, spike.preMaxCurve.t ./ 1000, spike.preMaxCurve.v, 'yo', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'y' )
  plot( ax, spike.postMaxCurve.t ./ 1000, spike.postMaxCurve.v, 'ys', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'y' )
  plot( ax, spike.maxDeriv.t ./ 1000, spike.maxDeriv.v, 'co', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'c' )
  plot( ax, spike.minDeriv.t ./ 1000, spike.minDeriv.v, 'cs', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'c')
  n1 = spike.n1List; n2 = spike.n2List;
  plot( ax, t(n1), v(n1), 'mo', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'm' )
  plot( ax, t(n2), v(n2), 'ms', ...
       'MarkerSize', markerSize, 'MarkerFaceColor', 'm' )
    
  plot( ax, spike.maxV.t ./ 1000, spike.maxV.v, 'ro', ...
        'MarkerSize', markerSize, 'MarkerFaceColor', 'r' )
  
  legendEntries = { 'pre max K', 'post max K', 'max dV', 'min dV', ...
                    'bracketStart', 'bracketEnd', 'maxV' };
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% add title and labels to axes
function addLabels( ax, titleStr, options )
  xlabel( ax, 'Time (s)', 'FontSize', options.labelFontSize )
  ylabel( ax, 'Voltage (mV)', 'FontSize', options.labelFontSize )
  title( ax, RealUnderscores( titleStr ), ...
         'FontSize', options.labelFontSize )
end