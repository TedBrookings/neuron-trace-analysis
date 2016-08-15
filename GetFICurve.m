% fI = GetFICurve(ephysData, options)
% Use GetSpikes.m to find the frequency and current injection for traces
%
% INPUTS:
%  -ephysData =  Data format from LoadAbf.m (see Get_traces.m)
% OPTIONS:
%  -vName        ['']    Name of voltage trace from ephysData. Autodetected
%                        if there is only one
%  -iName        ['']    Name of current trace. Autodetected if only one.
%  -plotSubject  [false] Don't plot if set to false or empty string. Make
%                        plots otherwise. If a string is provided, it will
%                        be incorporated into title/figure name
%  ADDITIONAL options are passed to GetSpikes.m, so run help(GetSpikes) to
%  see those options
% OUTPUTS:
%  -fI: a structure with fields giving information about the f-I curve
%    -f:  list of average spike frequency for each injection trace
%    -I:  list of mean current for each injection trace
%    -fLow: list of lowest frequency in confidence interval around f
%    -fHigh: list of highst frequency in confidence interval around f
%    -stdErrI: list of standard error of I
%    -IUnits: string specifying units of I ('nA' or 'pA')
%    -fitInfo: a structure specifying information obtained by fitting a
%              simple functinonal form to f-I (currently quadratic in I)
%       -IIntercept:  maximum current to yield f = 0
%       -fISlopeAtIIntercept: slope of the f-I curve at the I intercept
%       -fitCoefficients: coefficients obtained from fit
%       -fEstimator: a function that estimates f(I) when passed
%    -spikes: a structure containing spike detection info (the output of
%             call to GetSpikes.m)
function fI = GetFICurve(ephysData, varargin)
  parser = inputParser();
  parser.KeepUnmatched = true; % allow extra fields
  parser.addParameter( 'vName', '' )
  parser.addParameter( 'iName', '' )
  parser.addParameter( 'plotSubject', false )
  parser.addParameter( 'debugPlots', false )
  % make the I in FI curve be difference from holding current?
  parser.addParameter( 'subtractIHolding', true )
  % get options structure with all needed and extra fields
  parser.parse( varargin{:} )
  options = [struct2cell( parser.Results ); struct2cell( parser.Unmatched )];
  fNames = [ fieldnames( parser.Results ); fieldnames( parser.Unmatched ) ];
  options = cell2struct( options, fNames );

  %{
  % set the default options
  defaultOptions = { ...
    'vName', '', ...
    'iName', '', ...
    'plotSubject', false, ...
    'debugPlots', false ...
  };
  % get the options overrides from varargin
  options = GetOptions(defaultOptions, varargin, true);
  %}
  
  % get the voltage and current trace that correspond to the FI data we want
  [voltageTraceName, currentTraceName] = getTraceNames(ephysData, options);
  voltageTrace = ephysData.data.(voltageTraceName);
  currentTrace = ephysData.data.(currentTraceName);
  currentUnits = ephysData.units.(currentTraceName);
  
  % find spikes in the voltage trace
  dT = ephysData.time(2) - ephysData.time(1);
  spikes = getSpikes( voltageTrace, dT, options );
  
  % compute the array of injected currents from currentTrace
  % (use found spikes to remove artifacts in the current trace)
  [I, stdErrI, injectInds, noInjectInds] ...
    = getI(currentTrace, dT, spikes, options);
  
  % compute the array of spike frequencies and their uncertainties
  [numT, numTraces] = size(voltageTrace);
  [f, fLow, fHigh] = getF(dT, injectInds, spikes, numT, numTraces, options);
  
  % estimate IHolding, vHolding
  vHolding = voltageTrace(noInjectInds,:);
  deltaVHolding = std( vHolding(:) );
  vHolding = median( vHolding(:) );
  IHolding = currentTrace(noInjectInds,:);
  deltaIHolding = std( IHolding(:) );
  IHolding = median( IHolding(:) );
  
  if options.subtractIHolding
    I = I - IHolding;
  end
  
  % fit a simple functional form to the f-I curve
  fitInfo = fitCurve(f, I);
  
  % organize results into a structure
  fI = struct( ...
    'f', f, ...
    'I', I, ...
    'fLow', fLow, ...
    'fHigh', fHigh, ...
    'stdErrI', stdErrI, ...
    'vHolding', vHolding, ...
    'deltaVHolding', deltaVHolding, ...
    'IHolding', IHolding, ...
    'deltaIHolding', deltaIHolding, ...
    'IUnits', currentUnits, ...
    'fitInfo', fitInfo, ...
    'spikes', spikes ...
    );
  
  % make f-I plot if requested, do nothing if no plot requested
  plotFI(fI, options)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find the names of the intracellular voltage trace and the intracellular
% current trace
function [vName, iName] = getTraceNames(ephysData, options)
  % get the names of all the valid intracellular voltage or current traces
  vNames = {};
  iNames = {};
  traceNames = fieldnames( ephysData.units );
  for n = 1:numel( traceNames );
    name_n = traceNames{n};
    unit_n = ephysData.units.(name_n);
    if strcmp( unit_n, 'mV' )
      % trace has units of mV, is an intracellular voltage trace
      vNames = [vNames, name_n]; %#ok<AGROW>
    elseif ismember( unit_n, {'nA', 'pA'} )
      % trace has units of nA or pA, is an intracellular current trace
      iNames = [iNames, name_n]; %#ok<AGROW>
    end
  end
  
  % check to make sure that we found at least one of each type of trace
  if isempty( vNames )
    error( 'Could not find any intracellular voltage trace.' )
  elseif isempty( iNames )
    error( 'Could not find any intracellular current trace.' )
  end
  
  
  % filter by requested name
  if ~isempty( options.vName )
    match = strcmp( vNames, options.vName );
    if ~any( match )
      error( '%s does not match any voltage trace:%s\n', ...
             options.vName, sprintf( ' %s', vNames{:} ) )
    end
    vNames = vNames(match);
  end
  if ~isempty( options.iName )
    match = strcmp( iNames, options.iName );
    if ~any( match )
      error( '%s does not match any current trace:%s\n', ...
             options.iName, sprintf(' %s', iNames{:}) )
    end
    iNames = iNames(match);
  end
  
  % make sure that there is only one of each type of trace
  if numel( vNames ) > 1
    if stringBeginsWith( vNames{1}, 'Temp' )
      vNames(1) = [];
    elseif stringBeginsWith( vNames{2}, 'Temp' )
      vNames(2) = [];
    else
      error( 'Too many voltage traces: %s%s\n  (try specifying vName)', ...
             vNames{1}, sprintf( ', %s', vNames{2:end} ) )
    end
  end
  if numel( iNames ) > 1
    error( 'Too many current traces: %s%s\n  (try specifying iName)', ...
           iNames{1}, sprintf( ', %s', iNames{2:end} ) );
  end
  
  vName = vNames{1};
  iName = iNames{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function bw = stringBeginsWith( str, prefix )
  bw = numel( str ) >= numel( prefix ) ...
     && strcmp( str(1:numel( prefix )), prefix );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spikes = getSpikes( voltageTrace, dT, options )
  vDetect = reshape(voltageTrace, [], 1); % make into a column

  % do spike detection on vDetect
  removeOptions = {'vName', 'iName'};

  spikes = GetSpikes(dT, vDetect, rmfield(options, removeOptions));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute the array of injected currents from currentTrace
% This algorithm is a bit more complicated than obviously necessary,
% because the current may have large transient artifacts 
% (use found spikes to remove artifacts in the current trace)
function [I, stdErrI, injectInds, noInjectInds] ...
  = getI(currentTrace, dT, spikes, options)
  % remove artifacts from current trace due to spikes
  traceSize = size( currentTrace );
  currentTrace = reshape( currentTrace, [], 1 ); % make a column
  currentTrace = RemoveSpikes( dT, currentTrace, 'spikes', spikes, ...
                               'plot', options.debugPlots);
  currentTrace = reshape( currentTrace, traceSize );
  
  % find the current step with the largest 90th percentile current
  %[~, maxStepInd] = max( median(currentTrace, 1) );
  [~, maxStepInd] = max( quantile( currentTrace, 0.9, 1 ) );
  
  % assume that at least 5% of the time, the current is active. So find the
  % 95th percentile of current and use it as a preliminary threshold
  sortStep = sort( currentTrace(:, maxStepInd) );
  threshold = sortStep(round( 0.95 * numel( sortStep ) ));
  % find first and last indices current goes above this threshold
  i1 = find( currentTrace(:,maxStepInd) > threshold, 1 );
  i2 = find( currentTrace(:,maxStepInd) > threshold, 1, 'last' );
  
  % now assume that at least 50% of the time, the current is not in an
  % artifact. So find the 25th and 75th percentiles between i1 and i2
  sortStep = sort( currentTrace(i1:i2, maxStepInd) );
  lowThresh = sortStep(round( 0.25 * numel( sortStep ) ));
  highThresh = sortStep(round( 0.75 * numel( sortStep ) ));
  
  % average over indices where current is between those two thresholds
  injectInds = currentTrace(:,maxStepInd) >= lowThresh ...
    & currentTrace(:,maxStepInd) <= highThresh;
  
  % compute I as the median value of each column, on injectInds
  I = median( currentTrace(injectInds,:), 1 )';
  
  % compute standard error of I
  stdErrI = sqrt( var( currentTrace(injectInds,:), 1 )' ...
                  ./ sum( injectInds ) );
  
  % compute the indices when current is at holding value
  i1 = i1 - round( 25/dT );
  i2 = i2 + round( 25/dT );
  noInjectInds = [1:i1, i2:traceSize(1)];
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute spiking frequency and 1-sigma bounds
function [f, fLow, fHigh] = getF(dT, injectInds, spikes, numT, ...
                                 numTraces, options)
  % compute start and stop times of injection, (for computing rates)
  traceLen = numT * dT;
  tLow = dT * (find( injectInds, 1 ) - 1);
  tHigh = dT * (find( injectInds, 1, 'last' ) - 1);
  injectTime = 0.001 * (tHigh - tLow);
  
  % for each injection interval, compute f and upper/lower confidence
  % intervals
  f = zeros( numTraces, 1 );
  fLow = zeros( numTraces, 1 );
  fHigh = zeros( numTraces, 1 );
  for n = 1:numTraces
    spikeTimes_n = spikes.times(spikes.times >= tLow ...
                                & spikes.times < tHigh);
    % determine number of spikes used in estimate and the time interval used
    
    if isempty( spikeTimes_n )
      numSpikes = 0;
      tInterval = injectTime;
      stdErr = 1 / injectTime;  % assume +- 1 spike error
    elseif numel( spikeTimes_n ) == 1
      numSpikes = 1;
      tInterval = injectTime;
      stdErr = 1 / injectTime;  % assume +- 1 spike error
    else
      numSpikes = numel( spikeTimes_n ) - 1.0;
      tInterval = 0.001 * (max( spikeTimes_n ) - min( spikeTimes_n ));
      % compute stderr of mean of spike rates
      fList = 1000.0 ./ diff( spikeTimes_n );
      stdErr = std( fList ) / sqrt( numel( fList ) );
    end
    
    % compute frequencies
    f(n) = numSpikes / tInterval;
    fLow(n) = max( f(n) - stdErr, 0.0 );
    fHigh(n) = f(n) + stdErr;
    
    % update tLow and tHigh;
    tLow = tLow + traceLen;
    tHigh = tHigh + traceLen;
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% fit a simple functional form to the f-I curve
function fitInfo = fitCurve(f, I)
  % only fit the points where f > 0
  IFit = I(f > 0);
  fFit = f(f > 0);
  
  % return this if there's an error
  fitInfo = struct( ...
    'IIntercept', nan, ...
    'fISlopeAtIntercept', nan, ...
    'fitCoefficients', nan, ...
    'fEstimator', nan ...
  );

  if numel( IFit ) < 3
    % not enough info to do this fit
    return
  end

  % fit a quadratic polynomial to the first 4 non-zero points, unless there
  % are problems, in which case, increase the number of points used
  numPoints = min( 4, numel( IFit ) );
  coefs = polyfit( fFit(1:numPoints), IFit(1:numPoints), 2 );
  
  if coefs(1) <= 0
    %numPoints = numPoints + 1; % maybe loop this?
    numPoints = numel(fFit);
    coefs = polyfit( fFit(1:numPoints), IFit(1:numPoints), 2 );
  end
  if coefs(1) <= 0
    warning('WAVEFORM:FIFit', 'Quadratic fit has the wrong concavity!')
  end
  
  % This is a derivation of fMinQuad, the minimum of the parabola
  % I = coefs[1] * f^2 + coefs[2] * f + coefs[3]
  % dI = coefs[1] * 2 * f * df + coefs[2] * df
  % dI/df = 0 = coefs[1] * 2 * f + coefs[2]
  % fMin = -0.5 * coefs[2] / coefs[1]
  
  % compute IIntercept, fISlopeAtIntercept, and construct estimator of f(I)
  fMinQuad = -0.5 * coefs(2) / coefs(1);
  if fMinQuad < 0
    IIntercept = coefs(3);
    fISlopeAtIntercept = 1.0 / coefs(2);
  else
    IIntercept = coefs(3) + fMinQuad * (coefs(2) + fMinQuad * coefs(1));
    fISlopeAtIntercept = Inf;
  end
  
  if IIntercept > max(I)
    warning( 'Computed I-Intercept (%g) is greater than the maximum injected current (%g).', ...
             IIntercept, max(I))
    return
  end
  
  % construct an estimator of f(I)
  fEstimator = @(I) (I > IIntercept) .* ...
    ( sqrt( coefs(2)^2 + 4 * coefs(1) * (I - coefs(3)) ) ...
      - coefs(2)) / (2 * coefs(1) );
  
  fitInfo = struct( ...
    'IIntercept', IIntercept, ...
    'fISlopeAtIntercept', fISlopeAtIntercept, ...
    'fitCoefficients', coefs, ...
    'fEstimator', fEstimator ...
  );
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotFI(fI, options)
  if ischar(options.plotSubject)
    if isempty(options.plotSubject)
      % no plot desired
      return
    end
    plotTitle = [options.plotSubject, ': '];
  else
    if ~options.plotSubject
      % no plot desired
      return
    end
    plotTitle = '';
  end
  plotTitle = [plotTitle, 'FI Curve'];
  
  titleSize = 22;
  labelSize = 22;
  numFitPoints = 100;
  
  fig = NamedFigure( plotTitle ); fig.WindowStyle = 'docked'; clf( fig )
  ax = subplot( 1,1,1, 'Parent', fig );
  
  if ~isnan( fI.fitInfo.IIntercept )
    minPosI = min( fI.I(fI.I > fI.fitInfo.IIntercept) );
    IFit = [min(fI.I), fI.fitInfo.IIntercept, ...
            linspace( fI.fitInfo.IIntercept, minPosI, numFitPoints ), ...
            linspace( minPosI, max( fI.I ), numFitPoints )];
    fFit = fI.fitInfo.fEstimator(IFit);
    plot( ax, IFit, fFit, 'r-', 'LineWidth', 2 )
  end
  hold( ax, 'on' )
  plot( ax, fI.I, fI.f, 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 4 )
  for n = 1:numel( fI.f )
    plot( ax, [fI.I(n), fI.I(n)], [fI.fLow(n), fI.fHigh(n)], 'k-' )
  end
  for n = 1:numel( fI.I )
    ILow = fI.I(n) - fI.stdErrI(n);
    IHigh = fI.I(n) + fI.stdErrI(n);
    plot( ax, [ILow, IHigh], [fI.f(n), fI.f(n)], 'k-' )
  end
  hold( ax, 'off' )
  axis( ax, 'tight' )
  title( ax, plotTitle, 'FontName', 'Arial', 'FontSize', titleSize )
  xlabel( ax, sprintf('Injected Current (%s)', fI.IUnits ), ...
         'FontName', 'Arial', 'FontSize', labelSize )
  ylabel( ax, 'Spike Frequency (Hz)', 'FontName', 'Arial', ...
         'FontSize', labelSize )
end
