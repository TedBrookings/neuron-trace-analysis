% spike = GetSpikes(dT, v, plotSubject)
% Analyzes a single voltage waveform, looking for spikes
%    and bursts, and calculating relevant frequencies.
%
%  INPUT PARAMETERS:
%   -dT is sample time in ms
%   -v is array of voltages in mV
%    OPTIONAL:
%     -plotSubject should be set to true[false] to produce[suppress]
%       plots of waveforms/analysis.  Alternatively, it can be set
%       to a string to aid it titling plots (e.g. 'Exp #71')
%       plotSubject defaults to false
%     -lowCutoff: defaults to automatically detected. The threshold for
%       negative derivatives that constitutes a potential spike
%     -highCutoff: defaults to automatically detected. The threshold for
%       positive derivatives that constitutes a potential spike
%     -bracketWidth: defaults to 15ms. A spike must have a large positive
%       derivative followed by large negative within this interval
%     -minCutoffDiff: defaults to 0.1 (set to 0.001 for minis). If
%       autodetection produces high and low cutoffs less than this
%       difference, conclude there are no spikes.
%     -minSpikeHeight: default to 0.0 mV. Minimum allowable spike height to
%       be considered a valid spike.
%     -minSpikeAspect: defaults to 0.5 mV/ms. Minimum allowable ratio of
%       spike height to spike width to be considered a spike
%     -pFalseSpike: defaults to 0.05. Estimated proability of finding a
%       spurious spike in the whole trace
%     -recursive: defaults to false. if spikes are found, remove them and
%       try to find spikes in the remaining data. Keep doing this until no
%       new spikes are found
%     -debugPlots: defaults to false. When true, make extra plots depicting
%       the spike-finding process
%
%  OUTPUT PARAMETERS:
%   -spike:  a structure with the following fields
%    -spike.times is a plain list of spike times (in ms)
%    -spike.height is a plain list of spike heights (in mV)
%    -spike.width is a plain list of spike width (in ms)
%    -spike.freq is overall spiking frequency (in Hz)
%    -spike.intervals is a list of interspike intervals (in ms)
%    -spike.frequencies is a list of instantaneous frequencies (in Hz)
%    Shape information structures (should be self-descriptive)
%    -spike.maxV, spike.maxDeriv, spike.minDeriv, spike.preMinV,
%     spike.postMinV, spike.preMaxCurve, spike.postMaxCurve
%           Each contains a list of times/voltage points, and if relevant
%           another quantity (such as K for curvatures)
%
%List structures usually will have a name.list element, as well as
%  name.mean, name.stdDev, name.variance, name.coefOfVar
%  (a few are just plain lists)
%If a feature is not detected, relevant frequencies are set to
%  zero, and relevant lists are empty
%
function spike = GetSpikes(dT, v, varargin)
  if nargin < 2
    help GetSpikes
    error('Invalid number of arguments.')
  end
  if numel( dT ) > 1
    % user passed in array of time, rather than dT
    if numel( dT ) ~= numel( v )
      error( 'Time and Voltage arrays have different length!' )
    end
    dT = (dT(end) - dT(1)) / (length(dT) - 1);
  end
  
  if size( v,1 ) > 1
    if size( v,2 ) > 1
      error( 'Voltage must be a single array, not a matrix' )
    else
      v = v';
    end
  end

  % set the default options
  defaultOptions = { ...
    'plotSubject', false, ...
    'lowCutoff', NaN, ...
    'highCutoff', NaN, ...
    'bracketWidth', 3.0, ...
    'minCutoffDiff', 0.1, ...
    'minSpikeHeight', 4.0, ...
    'minSpikeAspect', 0.0, ...
    'noiseCheckQuantile', 0.67, ...
    'pFalseSpike', 1.0e-3, ...
    'distributionCheckProb', 0.5, ...
    'recursive', false, ...
    'discountNegativeDeriv', false, ...
    'removeOutliers', true, ...
    'outlierFraction', 0.33, ...
    'findMinis', false, ...
    'debugPlots', false ...
  };
  % get the options overrides from varargin
  [options, modified] = GetOptions( defaultOptions, varargin, true );

  if options.findMinis
    % if finding minis, change a few of the options (if not set by user)
    miniOptions = struct( ...
      'bracketWidth', 100.0, ...
      'minCutoffDiff', 0.001, ...
      'minSpikeHeight', 0.0, ...
      'minSpikeAspect', 0.0, ...
      'pFalseSpike', 0.25, ...
      'discountNegativeDeriv', true, ...
      'recursive', true, ...
      'outlierFraction', 0, ...
      'noiseCheckQuantile', 0.55 ...
     );
    for fName = fieldnames( miniOptions )'
      if ~modified.(fName{1})
        options.(fName{1}) = miniOptions.(fName{1});
      end
    end
  end

  %First get the spike times
  spike = getSpikeTimesThreshold( dT, v, options );
  if options.recursive
    oldSpikeTimes = [];
    while numel( oldSpikeTimes ) < numel( spike.times )
      oldSpikeTimes = spike.times;
      spike = getSpikeTimesThreshold( dT, v, options, spike );
    end
  end

  callstack = dbstack;
  if needPlot( options, callstack )
    hSpikes = PlotSpikes( dT, v, spike, [], options );
    
    % link relevant time axis together
    if options.debugPlots
      aSpikes = get(hSpikes, 'CurrentAxes');
      derivsTitle = makeTitle('Derivatives', options);
      %aDerivs = get(findobj('name', derivsTitle),'CurrentAxes');
      aDerivs = findobj('Tag', derivsTitle)';
      aHandles = [aSpikes, aDerivs];
      linkaxes(aHandles, 'x');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Finds spikes by looking for points where derivative is large
% (positive) followed quickly by a large (negative) derivative.
function spike = getSpikeTimesThreshold( dT, v, options, oldSpike )
  if dT < .005
    warning( 'WAVEFORM:SmallDT', ...
             'Very small dT (%g). Note dT should be in ms.', dT )
  end
  if nargin < 4
    oldSpike = [];
  end
  
  % get the voltage derivatives and thresholds for spike detection
  [deriv, deriv2, lowCutoff, highCutoff] = ...
    getDerivsAndThresholds( dT, v, options, oldSpike );
  
  % Get a list of putative spikes, bracketed between n1 and n2
  maxIndDiff = round( options.bracketWidth / dT );
  [n1List, n2List] = bracketSpikes( v, deriv, maxIndDiff, ...
                                    lowCutoff, highCutoff );
  
  %  Get spike shape
  spike = GetSpikeShape( n1List, n2List, dT, v, deriv, deriv2, options );
    
  %  Make plots if requested
  if needPlot(options) && options.debugPlots
    plotGetSpikeTimes( dT, v, deriv, deriv2, lowCutoff, highCutoff, ...
                       options );
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the voltage derivatives and thresholds for spike detection
function [deriv, deriv2, lowCutoff, highCutoff] = ...
    getDerivsAndThresholds(dT, v, options, oldSpike)
  maxTimeWidth = options.bracketWidth;
  nyquistRate = 1.0 / (2 * dT);
  fStop = min(nyquistRate * 2/3, 1.0 / maxTimeWidth);
  fPass = fStop;
  nyquistFrac = fStop / nyquistRate;
  [deriv, deriv2] = DerivFilter(v, dT, fPass, fStop);
  
  if isnan(options.lowCutoff) || isnan(options.highCutoff)
    [lowCutoff, highCutoff] = ...
      getAutoCutoffs(dT, deriv, nyquistFrac, options, oldSpike);
    if highCutoff - lowCutoff < options.minCutoffDiff
      % cutoffs are too closely spaced, corresponding to trivial spikes,
      % so widen them:
      fact = options.minCutoffDiff / (highCutoff - lowCutoff);
      highCutoff = highCutoff * fact;
      lowCutoff = lowCutoff * fact;
    end
  else
    if ~isnan(options.lowCutoff)
      lowCutoff = options.lowCutoff;
    end
    if ~isnan(options.highCutoff)
      highCutoff = options.highCutoff;
    end
  end
  
  if options.debugPlots
    titleStr = makeTitle('Spike Thresholds', options);
    fig = NamedFigure(titleStr); fig.WindowStyle = 'docked'; clf(fig)
    ax = subplot(1,2,1, 'Parent', fig);
    xRangeFull = [ min( deriv ), max( deriv ) ];
    xRange = [ max( 3 * lowCutoff, xRangeFull(1) ), ...
               min( 3 * highCutoff, xRangeFull(2) ) ];
    numInRange = sum( deriv >= xRange(1) & deriv <= xRange(2) );
    numBins = max(100, round( sqrt( numInRange) ));
    binW = diff( xRange ) / numBins;
    binMids = xRangeFull(1):binW:xRangeFull(2);
    [n, x] = hist(deriv, binMids);
    n = n ./ max(n);
    bar(ax, x, n, 1.0, 'EdgeColor', 'b', 'FaceColor', 'b');
    hold( ax, 'on' )
    plot(ax, [lowCutoff, lowCutoff], [0, 1], 'r')
    plot(ax, [highCutoff, highCutoff], [0, 1], 'g')
    hold(ax, 'off')
    xlabel(ax, 'Derivative (mV/ms)')
    ylabel(ax, 'Relative Frequency')
    title(ax, RealUnderscores(titleStr))
    legend( ax, { 'Derivatives', 'Low threshold', 'High threshold' }, ...
            'Location', 'Best' )
    axis( ax, 'tight' )
    xlim( ax, xRange )
    % we're debugging, so spit out information about the cutoffs
    fprintf('GetSpikes.m: low/high cutoff: %g/%g, bracketWidth=%g\n', ...
      lowCutoff, highCutoff, maxTimeWidth)
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get cutoffs for significant spiking
function [lowCutoff, highCutoff] = getAutoCutoffs(dT, deriv, ...
                                            nyquistFrac, options, oldSpike)

  if ~isempty( oldSpike )
    % first remove detected spikes from the list of voltage derivatives, then
    % sort into increasing order
    for n = numel( oldSpike.n1List ):-1:1
      n1 = oldSpike.n1List(n);
      n2 = oldSpike.n2List(n);
      deriv(n1:n2) = [];
    end
  end
  % sort the voltage derivative into a list of increasing order
  sortDeriv = sort( deriv(isfinite( deriv(:) )) );
  
  % number of *effective* trace points in a bracketed spike
  nBracket = nyquistFrac * options.bracketWidth / dT;
  % length of trace
  len = length(sortDeriv);
  logOdds = 4 * log(1 - options.pFalseSpike) / len / nBracket;
  
  % this is how rare a derivative has to be (either positive or negative) to
  % achieve the given false-detection probability
  minRareness = sqrt(-logOdds);
  
  % compute approximate 1/2-sigma levels for positive and negative
  % derivatives, based on presumably nearly-gaussian small derivatives near
  % the median derivative
  [peak, sigmaMinus, sigmaPlus] = FindPeak( sortDeriv, options.noiseCheckQuantile );
  wantedNumSigma = sqrt(2) * erfcinv( minRareness );
  highCutoff = max( [0, peak + wantedNumSigma * sigmaPlus] );
  
  if options.discountNegativeDeriv
    wantedNumSigma = min( 1, wantedNumSigma );
  end
  lowCutoff = min( [0, peak - wantedNumSigma * sigmaMinus] );
  %{ 
  highDV = sortDeriv(sortDeriv >= peak) - peak;
  highCutoff = peak + findThresh( highDV, minRareness, options );
  highCutoff = max(0, highCutoff);
  
  lowDV = flip( peak - sortDeriv(sortDeriv <= peak) );
  [lowThresh, lowSigma] = findThresh( lowDV, minRareness );
  if options.discountNegativeDeriv
    lowCutoff = peak - min(lowThresh, lowSigma);
  else
    lowCutoff = peak - lowThresh;
  end
  lowCutoff = min(0, lowCutoff);
  %}
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get 1-sided threshold for given rareness
function [thresh, sigma] = findThresh(data, rareness, options)
  checkP = options.noiseCheckQuantile;
  numData = numel(data);
  checkInd = 1 + round( (numData-1) * checkP );
  checkVal = data(checkInd);
  numSigmaCheck = sqrt(2) * erfcinv( checkP );
  sigma = checkVal / numSigmaCheck;
  
  wantedNumSigma = sqrt(2) * erfcinv( rareness );
  thresh = sigma * wantedNumSigma;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% high-pass filter a signal
function y = highPassFilter( y, halfFilterLength )
  % 1. prepare high-pass filter
  filtLen = 1 + 2 * halfFilterLength;
  filt = repmat( -1.0 / filtLen, 1, filtLen );
  filt(1 + halfFilterLength) = filt(1 + halfFilterLength) + 1.0;
  % 2. pad signal symmetrically
  y = [flip( y(2:halfFilterLength+1) ), ...
       y, ...
       flip( y(end-halfFilterLength-1:end-1) )];
  % 3. return valid part of convolution between padded-signal and filter
  y = conv( y, filt, 'valid' );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get a list of putative spikes, bracketed between n1 and n2
function [n1List, n2List] = bracketSpikes( v, deriv, maxIndDiff, ...
                                           lowCutoff, highCutoff )
% start looking for spikes at first sample where the derivative isn't very
% high
n1 = find(deriv < highCutoff, 1);
n1Barrier = 1;  % don't extend brackets past this number
numV = length(v);
n1Stop = numV - maxIndDiff;  % don't look past this barrier
n1List = [];
n2List = [];
while n1 < n1Stop
  if deriv(n1) < highCutoff
    n1 = n1 + 1;
  else  %Found potential beginning of a spike, try to bracket a spike
    n2 = n1 + 1;
    bracketSuccess = false;
    n2Stop = n1 + maxIndDiff;
    while n2 <= n2Stop
      if deriv(n2) > lowCutoff
        if deriv(n2) >= highCutoff
          % Slope is still high, reset n1
          n2Stop = min(n2, n1Stop) + maxIndDiff;
        end
        n2 = n2 + 1;
      else
        bracketSuccess = true;
        break
      end
    end
    if ~bracketSuccess
      n1 = n2 + 1;
      continue;
    end

    if n2 == numV || deriv(n2 + 1) > highCutoff || n2 - n1 < 2
      %probably just spurious
      n1 = n2 + 1;
      continue
    end

    %We've bracketed a spike between n1 and n2
    
    %We want to get some spike shape info, so extend n1 and n2
    %until we cross deriv = 0
    while n1 > n1Barrier && deriv(n1) > highCutoff
      n1 = n1 - 1;
    end
    while n2 < numV && deriv(n2) < lowCutoff
      n2 = n2 + 1;
    end
    
    n1List = [n1List, n1]; %#ok<AGROW>
    n2List = [n2List, n2]; %#ok<AGROW>
    n1Barrier = n2 + 1;
    n1 = n1Barrier;    
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function plotVar = needPlot(options, callStack)
if ischar(options.plotSubject)
  plotVar = true;
else
  plotVar = options.plotSubject;
end

if plotVar && nargin == 2 && length(callStack) >= 2
  plotVar = ~strcmp(callStack(2).name, 'AnalyzeWaveform');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot the derivatives and thresholds, showing how the affect spike
% detection
function fig = plotGetSpikeTimes( dT, v, deriv, deriv2, ...
                                  lowCutoff, highCutoff, options )
  titleStr = makeTitle('Derivatives', options);

  fig = NamedFigure( titleStr ); fig.WindowStyle = 'docked'; clf( fig )
  ax1 = subplot( 2, 1, 1, 'Parent', fig );
  
  numV = numel( v );
  dTSeconds = 0.001 * dT;
  tFinal = dTSeconds * (numV - 1);
  plot( ax1, 0:dTSeconds:tFinal, deriv, 'b-' )
  hold( ax1, 'on' )
  plot( ax1, [0, tFinal], [lowCutoff, lowCutoff], 'r-' )
  plot( ax1, [0, tFinal], [highCutoff, highCutoff], 'g-' )
  %xlabel( ax, 'Time (s)', 'FontSize', 1 8)
  ylabel( ax1, 'dV/dT (mV/ms)', 'FontSize', 18 )
  %title( ax, RealUnderscores( titleStr ), 'FontSize', 18 )
  legend( ax1, {'dV/dT', 'low threshold', 'high threshold'}, ...
         'Location', 'NorthOutside', 'Orientation', 'Horizontal' )
  hold( ax1, 'off' ) ; axis( ax1, 'tight' )
  ax1.Tag = titleStr;
  
  ax2 = subplot( 2, 1, 2, 'Parent', fig );
  plot( ax2, 0:dTSeconds:tFinal, deriv2, 'b-' )
  xlabel( ax2, 'Time (s)', 'FontSize', 18)
  ylabel( ax2, 'd^2V/dT^2 (mV/ms^2)', 'FontSize', 18 )
  axis( ax2, 'tight' )
  ax2.Tag = titleStr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set the full title for a figure based on base title and plotSubject
function titleStr = makeTitle( titleBase, options )
  if ischar(options.plotSubject)
    titleStr = [options.plotSubject, ': ', titleBase];
  else
    titleStr = titleBase;
  end
end
