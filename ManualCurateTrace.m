function waveInfo = ManualCurateTrace( traceInfo, spikeInfo, burstInfo, ...
                                        varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;
  parser.addParameter( 'dT', [] )
  parser.addParameter( 'title', 'Manually Curating Trace' )
  parser.addParameter( 'debugPlots', false )
  parser.addParameter( 'axesHeight', 0.85 )
  parser.addParameter( 'buttonWidth', 60 )
  parser.addParameter( 'buttonHeight', 40 )
  parser.addParameter( 'xPad', 5 )
  parser.addParameter( 'yPad', 5 )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  unmatchedFields = fieldnames( parser.Unmatched );
  for n = 1:numel( unmatchedFields )
    options.(unmatchedFields{n}) = parser.Unmatched.(unmatchedFields{n});
  end
  
  [t, v] = getTraceInfo( traceInfo, options );
  
  if ~exist( 'spikeInfo', 'var') || isempty( spikeInfo )
    spikeInfo = GetSpikes( t, v, 'plotSubject', options.debugPlots, ...
                           'debugPlots', options.debugPlots );
  end
  % spikeTimes: nSpikes x 3 matrix of times: spike start, peak, and end
  % spikePeaks: v(spike peak)
  [spikeTimes, spikePeaks] = getSpikeTimes( t, v, spikeInfo );

  if ~exist( 'burstInfo', 'var') || isempty( burstInfo )
    %burstInfo = GetBursts( t, v, spikeInfo );
    burstInfo = FindSpikeBursts( diff( t(1:2) ), v, spikeInfo );
  end
  % nBursts x 2 matrix of times: burst start, and end
  burstTimes = getBurstTimes( burstInfo );
  
  waveInfo = curateTrace( t, v, spikeTimes, spikePeaks, burstTimes, options );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, v] = getTraceInfo( traceInfo, options )
  if ischar( traceInfo )
    traceInfo = LoadAbf( traceInfo );
  end
  
  if isstruct( traceInfo )
    t = traceInfo.time; v = [];
    fNames = fieldnames( traceInfo.units );
    for n = 1:numel( fNames )
      if strcmp( traceInfo.units.(fNames{n}), 'mV' )
        v = traceInfo.data.(fNames{n});
        break
      end
    end
    if isempty( v )
      error( 'Could not find voltage trace in supplied data struct' )
    end
  elseif isfloat( traceInfo )
    if isempty( options.dT )
      error( 'When passing a trace, must supply dT as a parameter' )
    end
    v = traceInfo;
    t = options.dT .* (0:numel(v)-1);
  end
  
  if size( v, 2 ) ~= numel( t )
    v = v';
  end
  
  if ~isrow( v ) % v is a matrix
    [numTraces, numT] = size( v );
    dT = t(2) - t(1);
    % add a NaN to end of each epoch to draw a gap
    v = [v, NaN( numTraces, 1 ) ];
    tTemp = [t, t(end) + dT ];
    numT = numT + 1; epochDT = tTemp(end) + 0.25 * numT * dT;
    t(1, numel(v)) = 0;
    i2 = 0;
    for n = 1:numTraces
      i1 = i2 + 1; i2 = i2 + numT;
      t(i1:i2) = tTemp;
      tTemp = tTemp + epochDT;
    end
    v = v';
    v = v(:);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nSpikesx3 matrix of times: spike start, peak, and end
% spikePeaks: v(spike peak)
function [spikeTimes, spikePeaks] = getSpikeTimes( t, v, spikeInfo )
  spikeTimes = ...
    [t(spikeInfo.n1List)', t(spikeInfo.maxV.ind)', t(spikeInfo.n2List)'];
  spikePeaks = v(spikeInfo.maxV.ind);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nBursts x 2 matrix of times: burst start, and end
function burstTimes = getBurstTimes( burstInfo )
  burstTimes = [ burstInfo.startTime', burstInfo.stopTime' ];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create GUI for curation of traces
function traceInfo = curateTrace( t, v, spikeTimes, spikePeaks, ...
                                  burstTimes, options )
  % create GUI figure
  titleStr = options.title;
  fig = NamedFigure( titleStr ); clf( fig ) ; fig.WindowStyle = 'docked';

  % store trace data that is necessary for curation
  traceData = struct( 'dT', diff( t(1:2) ), 't', t, ...
                      'v', v, 'minV', min( v ), 'maxV', max( v ),...
                      'spikeTimes', spikeTimes, 'spikePeaks', spikePeaks,...
                      'burstTimes', burstTimes );
  traceData.originalTraceData = traceData;
  fig.UserData = traceData;

  % add axes to figure
  addAxes( fig, options );
    
  % add control panel to figure
  addControlPanel( fig, ax, options )

  % add resize function
  fig.ResizeFcn = @(varargin) figResizeCallback( ax, options );

  % plot the trace
  plotTrace( ax )
  
  % wait for curation to finish
  uiwait( fig )

  % process the raw trace data
  fig.UserData = getProcessedTraceInfo( fig.UserData );

  % return final curated trace information
  traceInfo = fig.UserData;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add axes to figure
function addAxes( fig, options )
  h = options.axesHeight; y = 1.0 - h;
  ax = axes( 'Parent', fig, 'Units', 'Normalized', 'Position', [0 y 1 h], ...
             'Tag', 'Curate Trace Axes' );
  ax.Units = 'Points';
  yChars = max( cellfun( @(l) numel( l ), ax.YTickLabel ) );
  if isempty( yChars ) || yChars < 1, yChars = 1; end
  ax.Position = ax.Position + ax.FontSize .* [yChars, 2, -yChars - 1, -3];
  ax.Units = 'normalized';
  hold( ax, 'on' )
  ax.UserData = struct( 'status', 'normal', 'editX', NaN );
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% add control panel
function addControlPanel( fig, ax, options )
  buttonWidth = options.buttonWith; buttonHeight = options.buttonHeight;
  xPad = options.xPad; yPad = options.yPad;

  h = ax.Position(2) - yPad;
  controlPanel = uipanel( fig, 'Position', [ 0, 0, 1, h], ...
                          'Tag', 'ControlPanel' );
  % add okay button
  x = xPad; y = yPad;
  uicontrol( controlPanel, 'String', 'Okay', 'Callback', @okayCallback, ...
             'Style', 'pushbutton', 'Units', 'pixels', ...
             'Position', [x y buttonWidth buttonHeight], ...
             'BackgroundColor', [0.5 1 0.5]);
  % add reset button
  x = x + buttonWidth + 2 * xPad;
  uicontrol( controlPanel, 'String', 'Reset', 'Callback', @resetCallback, ...
             'Style', 'pushbutton', 'Units', 'pixels', ...
             'Position', [x y buttonWidth buttonHeight] );
  % add clear button
  x = x + buttonWidth + 2 * xPad;
  uicontrol( controlPanel, 'String', 'Clear', 'Callback', @clearCallback, ...
             'Style', 'pushbutton', 'Units', 'pixels', ...
             'Position', [x y buttonWidth buttonHeight] );

  % add instructions button
  x = x + buttonWidth + 2 * xPad;
  uicontrol( controlPanel, 'String', 'Instructions', ...
             'Callback', @instructionsCallback, ...
             'Style', 'pushbutton', 'Units', 'pixels', ...
             'Position', [x y buttonWidth buttonHeight] );
end
                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figResizeCallback( ax, options )
  h = options.axesHeight; y = 1.0 - h;
  ax.Position = [0 y 1 h];
  ax.Units = 'Points';
  yChars = max( cellfun( @(l) numel( l ), ax.YTickLabel ) );
  if isempty( yChars ) || yChars < 1, yChars = 1; end
  ax.Position = ax.Position + ax.FontSize .* [yChars, 2, -yChars - 1, -3];
  ax.Units = 'normalized';
  
  fig = ax.Parent;
  controlPanel = findobj( fig.Children, 'Tag', 'ControlPanel' );
  controlPanel.Position(4) = ax.Position(2) - options.yPad;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotTrace( ax )
  fig = ax.Parent; traceData = fig.UserData;
  
  cla( ax ); hold( ax, 'on' )
  plot( ax, traceData.t .* 1e-3, traceData.v )
  plot( ax, traceData.spikeTimes(:,2) .* 1e-3, traceData.spikePeaks, 'co' )
  h = traceData.maxV - traceData.minV;
  for n = 1:size( traceData.burstTimes, 1 )
    t1 = traceData.burstTimes(n,1) / 1000;
    t2 = traceData.burstTimes(n,2) / 1000;
    rectangle( 'Position', [t1, traceData.minV, t2 - t1, h], ...
               'FaceColor', [1 0 0 0.2], 'EdgeColor', [1 0 0] )
  end
  
  ax.HitTest = 'on';
  ax.PickableParts = 'all';
  for n = 1:numel( ax.Children )
    ax.Children(n).HitTest = 'off';
  end

  ax.ButtonDownFcn = @axClickCallback;
  ax.Tag = 'Curate Trace Axes';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% callback to execute when curation axis is clicked on
function axClickCallback( ax, eventData )
  % get the x,y position of the click
  clickX = eventData.IntersectionPoint(1);
  clickY = eventData.IntersectionPoint(2);
  % get the button number:
  %  -1 (left click / regular press on touch screen)
  %  -3 (right click / long press on touch screen)
  buttonNumber = eventData.Button;

  % get the selection type:
  %  -normal
  %  -extend (shift was pressed)
  %  -alt (Ctrl was pressed, or right-clicking)
  %  -open (double-click)
  fig = ax.Parent;
  selectionType = fig.SelectionType;

    % get the edit status
  %  -{'normal', nan} normal
  %  -{'spike', int} one side of spike bracket has been selected, and index
  %                  of that selection
  %  -{'burst', int} one side of burst has been selected, and index of that
  %                  selection
  editStatus = ax.UserData;

  traceData = fig.UserData;
  switch buttonNumber
    case 1 % left click / regular press on touch screen
      % add a spike or burst
      switch editStatus.status
        case 'normal'
          % no selection in progress, so decide what to do based on
          % selectionType
          switch selectionType
            case 'normal' % try to bracket spike near this one click
              % get the closest point on the trace to this click
              [minDist, closeInd] = findClosestPoint( ax, traceData, ...
                                                      clickX, clickY );
              if minDist > 0.1, return ; end
              %{
              garbage to be moved:
              if isempty( traceData.curated.originalSpikes )
                traceData.curated.originalSpikes = traceData.spikes;
              end
              %}
              traceData = addSpike( traceData, closeInd );
            case 'alt' % select one side of spike bracket
              % update editStatus
              ax.UserData = struct( 'status', 'spike', 'editX', clickX );
              return
            case 'extend' % select one side of burst
              % update editStatus
              ax.UserData = struct( 'status', 'burst', 'editX', clickX );
              return
            otherwise % double click, ignore
              return
          end
          
        case 'spike'
          % one side of spike bracket is selected, so use this click as the
          % other side
          traceData = addSpike( traceData, editStatus.editX, clickX );
          % set editStatus to normal
          ax.UserData = struct( 'status', 'normal', 'editX', NaN );
        case 'burst'
          % one side of burst is selected, so use this click as the other
          % side
          traceData = addBurst( traceData, editStatus.editX, clickX );
          % set editStatus to normal
          ax.UserData = struct( 'status', 'normal', 'editX', NaN );
        otherwise % this should be impossible
          error( 'Invalid edit status: %s', editStatus )
      end
      
    case 3 % right click / long press on touch screen
      % remove a spike or burst, or cancel current selection
      if strcmp( editStatus, 'normal' )
        % no selection in progress, delete based on selectionType
        switch selectionType
          case 'open' % ignore double-click
            return
          case 'extend' % delete nearest burst
            traceData = removeNearestBurst( ax, traceData, clickX );
          otherwise % delete nearest spike
            traceData = removeNearestSpike( ax, traceData, ...
                                            clickX, clickY );
        end
      else
        % spike or burst was selected, just cancel the selection
        % set editStatus to normal
        ax.UserData = struct( 'status', 'normal', 'editX', NaN );
        return
      end
      
    otherwise
      fprintf( 2, 'Don''t press button %d\n', eventData.button );
      return
  end

  % update traceData
  fig.UserData = traceData;
  
  % update plot of trace
  plotTrace( ax )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the closest point on vTrace to this click, and the distance
function [minDist, closeInd] = findClosestPoint( axesObj, traceData, ...
                                                 clickX, clickY, t, v )
  xRange = 1.0e3 * diff( xlim( axesObj ) );yRange= diff( ylim( axesObj ) );
  
  % get square distance to click
  if ~exist( 't', 'var' ), t = traceData.t; end
  if iscolumn( t ), t = t'; end
  if ~exist( 'v', 'var' ), v = traceData.v; end
  if iscolumn( v ), v = v'; end
  clickDist = ( (t - clickX * 1.0e3) ./ xRange ).^2 ...
            + ( (v - clickY) ./ yRange ).^2;
  % find closest point
  [minDist, closeInd] = min( clickDist );
  minDist = sqrt( minDist );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function traceData = addBurst( traceData, t1, t2 )
  % ensure t1 < t2
  if t2 < t1, tTemp = t1; t1 = t2; t2 = tTemp; end

  % insert this burst into traceInfo.burstTimes
  burstInd = find( traceData.burstTimes(:,1) < t1, 1, 'last' );
  if isempty( burstInd ), burstInd = 0; end
  traceData.burstTimes = [ traceData.burstTimes(1:burstInd,:) ; ...
                           [t1, t2] ; ...
                           traceData.burstTimes(burstInd+1:end,:) ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function traceData = addSpike( traceData, t1, t2 )
  if nargin < 3 || isempty( t2 )
    % note, in this case t1 is actually an index to .t and .v
    [t1, t2] = bracketSpike( traceData, t1 );
  else
    % ensure t1 < t2
    if t2 < t1, tTemp = t1; t1 = t2; t2 = tTemp; end
  end
  
  % find time and voltage of spike peak
  i1 = find( traceData.t >= t1, 1 );
  i2 = i1 + find( traceData.t(i1+1:end) <= t2, 1, 'last' );
  [peakV, peakInd] = max( traceData.v(i1:i2), [], 'omitnan' );
  spikeT = traceData.t(peakInd);
  
  % insert this spike into traceData
  spikeInd = find( traceData.spikeTimes(:,2) < spikeT, 1, 'last' );
  if isempty( spikeInd ), spikeInd = 0; end
  traceData.spikeTimes = [ traceData.spikeTimes(1:spikeInd,:) ; ...
                           [t1, spikeT, t2] ; ...
                           traceData.spikeTimes(spikeInd+1:end,:) ];
  traceData.spikePeaks = [ traceData.spikePeaks(1:spikeInd) ; ...
                           peakV ; ...
                           traceData.spikePeaks(spikeInd+1:end) ];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t1, t2] = bracketSpike( traceData, indClick )
  % first find direction towards peak
  v = traceData.v; numV = numel( v );
  iRange = [ max( indClick - 10, 1 ), min( indClick + 10, numV ) ];
  [~, iPeak] = max( v(iRange), [], 'omitnan' );
  iPeak = iRange(iPeak);
  if iPeak < indClick
    i1 = max( 1, iPeak - 1 ); i2 = indClick;
  elseif iPeak > indClick
    i1 = indClick ; i2 = min( numV, iPeak + 1 );
  else
    i1 = max( 1, iPeak - 1 ); i2 = min( numV, iPeak + 1 );
  end
  
  % extend i1 to the left until spike is bracketed
  while i1 > 1
    if v(i1-1) < v(i1)
      i1 = i1 - 1;
    elseif i1 > 2 && v(i1 - 2) < v(i1)
      i1 = i1 - 2;
    else
      break
    end
  end
  
  % extend i2 to the right until spike is bracketed
  while i2 < numV
    if v(i2+1) < v(i2)
      i2 = i2 + 1;
    elseif i2 < numV - 1 && v(i2 +2) < v(i2)
      i2 = i2 + 2;
    else
      break
    end
  end
  
  % convert indices to times
  t1 = 1.0e-3 * traceData.t(i1);
  t2 = 1.0e-3 * traceData.t(i2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove spike nearest to click
function traceData = removeNearestSpike( ax, traceData, clickX, clickY)
  % get the closest spike peak to this click  
  t = traceData.spikeTimes(:,2);
  if isempty( t ), return, end
  v = traceData.spikePeaks;
  [minDist, removeInd] = findClosestPoint( ax, traceData, clickX, clickY,...
                                           t, v );
  % if click is far away from any spike, return without change
  if minDist > 0.1, return , end
  % remvoe selected spike
  traceData.spikeTimes(removeInd,:) = [];
  traceData.spikePeaks(removeInd) = [];  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove burst nearest to click
function traceData = removeNearestBurst( traceData, clickX )
  t = traceData.burstTimes;
  if isempty( t ), return, end
  
  clickX = clickX * 1.0e3;
  rmInd = find( t(:,1) <= clickX & clickX <= t(:,2), 1 );
  if ~isempty( rmInd )
    traceData.burstTimes(:,rmInd) = [];
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function okayCallback( buttonObj, eventData )
  uiresume()
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function resetCallback( buttonObj, eventData )
  % restore original traceData;
  fig = buttonObj.Parent.Parent;
  traceData = fig.UserData;
  originalTraceData = traceData.originalTraceData;
  traceData = rmfield( traceData, 'originalTraceData' );
  if isequaln( traceData, originalTraceData )
    return % nothing to do
  end
  originalTraceData.originalTraceData = originalTraceData;
  fig.UserData = originalTraceData;

  % update plot of trace
  ax = findobj( fig.Children, 'Tag', 'Curate Trace Axes' );
  plotTrace( ax )    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function clearCallback( buttonObj, eventData )
  % remove burst times
  fig = buttonObj.Parent.Parent;
  if isempty( fig.UserData.burstTimes )
    return % no need to do anything
  end
  fig.UserData.burstTimes = [];
  
  % update plot of trace
  ax = findobj( fig.Children, 'Tag', 'Curate Trace Axes' );
  plotTrace( ax )  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function instructionsCallback( buttonObj, eventData )
  instructions = {...
    'APPROVING SOURCES:', ...
    '  Click "Okay" button to approve spikes and bursts', ...
    '  Click "Clear" button to clear all bursts', ...
    '  Click "Reset" button to reset spikes and bursts to auto-detected', ...
    '    (only available if source is already approved/rejected)', ...
    '  Click "Quit" button to save and quit', ...
    'ANNOTATING SOURCES:', ...
    '  Left click near trace to add a spike.', ...
    '  To forceably bracket spike, CTRL-left-click in two locations:', ...
    '    -spike will be in between those locations', ...
    '    -to abort adding spike, right click', ...
    '  To add a burst, SHIFT-left-click in two locations:', ...
    '    -burst will be in between those locations', ...
    '    -to abort adding burst, right click', ...
    '  To delete a spike, right click', ...
    '  To delete a burst, SHIFT-right-click' ...
    };
  msgbox( instructions )
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% process the raw trace data
function traceInfo = getProcessedTraceInfo( traceData )
  dT = traceData.dT; v = traceData.v; t = traceData.t;
  n1List = 1 + round( (traceData.spikeTimes(:,1)' - t(1))  ./ dT );
  n2List = 1 + round( (traceData.spikeTimes(:,3)' - t(1))  ./ dT );
  spike = GetSpikeShape( n1List, n2List, dT, v, [], [], ...
                         'removeOutliers', false, 'pFalseSpike', 1.99 );
  
                       
  traceInfo = struct( ...
    'dT', dT, 't', t, 'v', v, ...
    'spike', spike, 'burst', burst );
end