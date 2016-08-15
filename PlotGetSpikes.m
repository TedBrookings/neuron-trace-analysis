function varargout = PlotGetSpikes(dt, v, spike, options, burst)
% PlotGetSpikes(dt, v, spike, options, burst)
% plot spike info (and optionally bursts) on top of voltage trace
%  INPUTS:
%   -t:  time array (ms)
%   -v:  voltage trace array (mV)
%   -spike: structure with spike information, as calculated by GetSpikes.m
%   -options: structure with detection/plot options
%   OPTIONAL:
%   -burst:  structure with burst information, as caculated by AnalyzeBurst.m

if length(dt) == 1
  % passed dt
  dt = dt / 1000; % convert to sec
  tFinal = dt * (length(v) - 1);
  t = 0:dt:tFinal;
else
  % passed t (as dt)
  t = 0.001 * dt;
  if t(1) ~= 0
    t = t - t(1);
  end
  dt = (t(end) - t(1)) / (length(v) - 1);
end
  
spikeTimes = spike.times / 1000;
if nargin == 5
  % using burst information
  burstTimes = burst.Times / 1000;
  burstLen = burst.Durations.list / 1000;
  baseTitle = 'Spike/Burst Detection';
else
  % only using spike information
  baseTitle = 'Spike Detection';
end

% define upper and lower voltages of lines indicating spike times
top = max(v);
bottom = min(v);
delta = 0.1 * (top - bottom);
bottom = bottom - delta;
top = top + delta;

% make the title
if ischar(options.plotSubject)
  titleStr = [options.plotSubject, ': ', baseTitle];
else
  titleStr = baseTitle;
end

if ~isfield(options, 'timesOnly')
  options.timesOnly = false;
end
if ~options.timesOnly
  % get some shape information
  if isfield(spike.maxV.v, 'list')
    % spike is "structified"
    vMax = spike.maxV.v.list;
    vPreMaxK = spike.preMaxCurve.v.list;
    vPostMaxK = spike.postMaxCurve.v.list;
    vPreMinV = spike.preMinV.v.list;
    vPostMinV = spike.postMinV.v.list;
    vMaxDeriv = spike.maxDeriv.v.list;
    vMinDeriv = spike.minDeriv.v.list;
  else
    vMax = spike.maxV.v;
    vPreMaxK = spike.preMaxCurve.v;
    vPostMaxK = spike.postMaxCurve.v;
    vPreMinV = spike.preMinV.v;
    vPostMinV = spike.postMinV.v;
    vMaxDeriv = spike.maxDeriv.v;
    vMinDeriv = spike.minDeriv.v;
  end
end 

% create a new figure, name it, and make it ready for plotting
h = NamedFigure(titleStr);
set(h, 'WindowStyle', 'docked');
clf; % clear the plot, in case something is already there
whitebg(h, 'k');  % make the background black
hold on;

legendEntries = {};
if nargin == 5
  % first draw blue rectangle to signify burst times
  numBursts = length(burstTimes);
  if numBursts > 0
    burstGroup = hggroup;
    for n = 1:numBursts
      tLow = burstTimes(n);
      tHigh = tLow + burstLen(n);
      fill([tLow, tHigh, tHigh, tLow], [bottom, bottom, top, top], 'b', ...
        'Parent', burstGroup);
    end
    set(get(get(burstGroup,'Annotation'),'LegendInformation'),...
      'IconDisplayStyle','on');
    legendEntries = [legendEntries, 'Burst'];
  end
end

numSpikes = length(spikeTimes);

if options.timesOnly
  % no shape information, just mark spikes
  for n=1:numSpikes
    % overlay red lines indicating spikes
    plot([spikeTimes(n), spikeTimes(n)], [bottom,top], 'r-');
  end
  
  % finally draw the voltage trace in white:
  plot(t, v, 'w-');
else
  % we have shape information, show some of it
  
  % draw the voltage trace in white:
  voltageLine = plot(t, v, 'w-');
  set(get(get(voltageLine,'Annotation'),'LegendInformation'),...
    'IconDisplayStyle','off'); % Exclude line from legend

  plot( spike.preMaxCurve.t ./ 1000, vPreMaxK, 'yo', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'y' )
  plot( spike.postMaxCurve.t ./ 1000, vPostMaxK, 'ys', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'y' )
  plot( spike.maxDeriv.t ./ 1000, vMaxDeriv, 'co', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'c' )
  plot( spike.minDeriv.t ./ 1000, vMinDeriv, 'cs', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'c')
  n1 = spike.n1List; n2 = spike.n2List;
  plot( t(n1), v(n1), 'mo', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'm' )
  plot( t(n2), v(n2), 'ms', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'm' )
  
  plot( spike.maxV.t ./ 1000, vMax, 'ro', ...
        'MarkerSize', 6, 'MarkerFaceColor', 'r' )
end


if numSpikes > 0
  if options.timesOnly
    legendEntries = [legendEntries, 'spike'];
  else
    legendEntries = [legendEntries, ...
      { 'pre max K', 'post max K', 'max dV', 'min dV', ...
        'bracketStart', 'bracketEnd', 'maxV' } ];
  end
end
if ~isempty(legendEntries)
  legend(legendEntries{:}, 'Location', 'SouthOutside', 'Orientation', 'Horizontal')
end

% axis labels and title
xlabel('Time (s)', 'FontSize', 18)
ylabel('Voltage (mV)', 'FontSize', 18)
title(RealUnderscores(titleStr), 'FontSize', 18)
hold off;
if nargout
  varargout = {h};
end
return
