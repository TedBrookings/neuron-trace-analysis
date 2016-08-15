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
  if length(dT) > 1
    % user passed in array of time, rather than dT
    if length(dT) ~= length(v)
      error('Time and Voltage arrays have different length!')
    end
    dT = (dT(end) - dT(1)) / (length(dT) - 1);
  end
  
  if size(v,1) > 1
    if size(v,2) > 1
      error('Voltage must be a single array, not a matrix')
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
    'minSpikeHeight', 0.0, ...
    'minSpikeAspect', 0.0, ...
    'pFalseSpike', 1.0e-2, ...
    'distributionCheckProb', 0.5, ...
    'recursive', false, ...
    'discountNegativeDeriv', false, ...
    'removeOutliers', true, ...
    'findMinis', false, ...
    'debugPlots', false ...
  };
  % get the options overrides from varargin
  [options, modified] = GetOptions(defaultOptions, varargin, true);

  if options.findMinis
    % if finding minis, change a few of the options (if not set by user)
    miniOptions = struct( ...
      'bracketWidth', 50.0, ...
      'minCutoffDiff', 0.001, ...
      'minSpikeAspect', 0.0, ...
      'pFalseSpike', 0.05, ...
      'discountNegativeDeriv', true, ...
      'recursive', true ...
     );
    for fName = fieldnames(miniOptions)'
      if ~modified.(fName{1})
        options.(fName{1}) = miniOptions.(fName{1});
      end
    end
  end

  %First get the spike times
  spike = getSpikeTimesThreshold(dT, v, options);
  if options.recursive
    oldSpikeTimes = [];
    while length(oldSpikeTimes) < length(spike.times)
      oldSpikeTimes = spike.times;
      spike = getSpikeTimesThreshold(dT, v, options, spike);
    end
  end

  %Next get the overall spike frequency
  spike.freq = getSpikeFrequency(spike.times, dT * (length(v) - 1));

  callstack = dbstack;
  if needPlot(options, callstack)
    hSpikes = PlotGetSpikes(dT, v, spike, options);
    
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
function spike = getSpikeTimesThreshold(dT, v, options, oldSpike)
  if dT < .005
    warning('WAVEFORM:SmallDT', ...
            'Very small dT (%g). Note dT should be in ms.', dT)
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
  spike = getSpikeShape(n1List, n2List, dT, v, deriv, deriv2, options);
  
  %  Calculate spike intervals and frequencies
  if isempty(spike.times)
    spike.intervals = [];
    spike.frequencies = [];
  else
    spike.intervals = spike.times(2:end) - spike.times(1:(end-1));
    spike.frequencies = 1000 ./ spike.intervals;
  end
  
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
    numBins = max(100, round( sqrt( numel(deriv) ) ));
    [n, x] = hist(deriv, numBins);
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
    xRange = xlim();
    xRange = [ max( 3 * lowCutoff, xRange(1) ), ...
               min( 3 * highCutoff, xRange(2) ) ];
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

  if ~isempty(oldSpike)
    % first remove detected spikes from the list of voltage derivatives, then
    % sort into increasing order
    for n = length(oldSpike.n1List):-1:1
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
  peak = findPeak( sortDeriv );
  
  highDV = sortDeriv(sortDeriv >= peak) - peak;
  highCutoff = peak + findThresh( highDV, minRareness );
  highCutoff = max(0, highCutoff);
  
  lowDV = flip( peak - sortDeriv(sortDeriv <= peak) );
  [lowThresh, lowSigma] = findThresh( lowDV, minRareness );
  if options.discountNegativeDeriv
    lowCutoff = peak - min(lowThresh, lowSigma);
  else
    lowCutoff = peak - lowThresh;
  end
  lowCutoff = min(0, lowCutoff);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get 1-sided threshold for given rareness
function [thresh, sigma] = findThresh(data, rareness)
  checkP = 0.5;
  numData = numel(data);
  checkInd = 1 + round( (numData-1) * checkP );
  checkVal = data(checkInd);
  numSigmaCheck = sqrt(2) * erfcinv( checkP );
  sigma = checkVal / numSigmaCheck;
  
  wantedNumSigma = sqrt(2) * erfcinv( rareness );
  thresh = sigma * wantedNumSigma;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quick method of estimating the height of noise in the trace, old method
function noiseHeight = getNoiseHeightFast(v, n1List, n2List, options)
  if isempty( n1List )
    noiseHeight = 0;
    return
  end
  % eliminate spikes from consideration;
  spikeInds = arrayfun( @(n1,n2) n1:n2, n1List, n2List, ...
                        'UniformOutput', false );
  spikeInds = cat(2, spikeInds{:});
  minNumDataPoints = 100; %need this many data points to make an okay guess
  if numel(spikeInds) > numel(v) - minNumDataPoints
    % not a lot of non-spike data to work with.
    % assume spikes should be more than just large single-point fluctuations,
    % so we want outliers on individual point-to-point differences
    noiseHeights = abs( diff( v ) );
    
    % find threshold such that an individual noiseHeight is highly unlikely
    % to exceed it in the whole movie
    noiseHeights = noiseHeights(noiseHeights > 0);
    noiseHeights = sort(noiseHeights);
    % find the peak of those noise heights
    peakNoise = findPeak( noiseHeights );
    highNoiseHeights = noiseHeights(noiseHeights > peakNoise) - peakNoise;
    % assume peak is ~ gaussian, and estimate sigma of that peak by finding the
    % location halfway down the cumulative distribution
    numSigmaHalf = sqrt(2.0) * erfinv(0.5);
    sigma = median(highNoiseHeights) / numSigmaHalf;
    % choose threshold so rare the the probability of a false spike in whole data
    % set is pFalseSpike. do numerically more stable version of this:
    % rareness = 1 - (1 - pFalseSpike).^(1.0 / numel(noiseHeights));
    rareness = -expm1(log1p(-options.pFalseSpike)) / numel(highNoiseHeights);
    numSigmaNeeded = sqrt(2) * erfcinv(rareness);
    
    noiseHeight = peakNoise + sigma * numSigmaNeeded;
  else
    % enough non-spike data. Try to estimate the noise as 2 * the standard
    % deviation of high-pass filtered non-spike data
    % high-pass filter v
    spikeWidth = round(median(n2List - n1List));
    halfLen = round( (spikeWidth - 1) / 2 );
    noiseHeights = highPassFilter(v, halfLen);
    noiseHeights( spikeInds ) = [];
    numSigmaNeeded = sqrt(2) * erfcinv( options.pFalseSpike );
    noiseHeight = numSigmaNeeded * std(noiseHeights);
  end
  
  if options.debugPlots
    titleStr = makeTitle('Spike Thresholds', options);
    fig = NamedFigure(titleStr); fig.WindowStyle = 'docked';
    ax = subplot(1,2,2, 'Parent', fig);
    numBins = max(100, round( sqrt( numel(noiseHeights) ) ));
    [n, x] = hist(noiseHeights, numBins);
    n = n ./ max(n);
    bar(ax, x, n, 1.0, 'EdgeColor', 'b', 'FaceColor', 'b');
    hold( ax, 'on' )
    plot(ax, [noiseHeight, noiseHeight], [0, 1], 'g')
    hold(ax, 'off')
    xlabel(ax, 'Noise (mV)')
    ylabel(ax, 'Relative Frequency')
    titleStr = makeTitle('Spike Height Threshold', options);
    title(ax, RealUnderscores(titleStr))
    legend(ax, 'Noise', 'Spike height threshold', 'Location', 'Best')
    % we're debugging, so spit out information about the cutoffs
    fprintf( 'GetSpikes.m: spike height cutoff: %g\n', ...
             noiseHeight )
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quick method of estimating the height of noise in the trace
function noiseHeight = getNoiseHeight(v, n1List, n2List, options)
  if isempty( n1List )
    noiseHeight = 0;
    return
  end
  
  % find indices when a spike is in progress
  spikeInds = arrayfun( @(n1,n2) n1:n2, n1List, n2List, ...
                        'UniformOutput', false );
  spikeInds = cat(2, spikeInds{:});
  % how wide are spikes
  spikeWidth = median( n2List - n1List );
  
  minNumDataPoints = 100; %need this many data points to make an okay guess
  fewDataPoints = numel( spikeInds ) > numel( v ) - minNumDataPoints;
  if fewDataPoints
    % not a lot of non-spike data to work with.
    % assume spikes should be more than just large single-point fluctuations,
    % so we want outliers on individual point-to-point differences
    noiseHeights = abs( diff( v ) );
  else % enough spike data
    % get noise heights as v - lower envelope (v )
    % set filter length to odd integer ~ 1/2 spike width
    filtLen = 1 + 2 * round( (spikeWidth - 1) / 4 );
    % get lower-envelope of trace
    envelope = LowerEnvelope( v, filtLen, 'plot', options.debugPlots, ...
                          'title', makeTitle( 'LowerEnvelope', options ) );
    % noise heights are the height above envelope
    noiseHeights = v - envelope;
    %remove spike indices from consideration
    noiseHeights(spikeInds) = [];
    
    zeroInds = find( noiseHeights == 0 );
    noiseHeights = arrayfun( @(i1,i2) max( noiseHeights(i1:i2) ), ...
                             [1, zeroInds], ...
                             [zeroInds, numel( noiseHeights )] );
  end

  noiseHeights = noiseHeights(noiseHeights(:) > 0);
  noiseHeights = sort(noiseHeights);
  % find the peak of those noise heights
  peakNoise = findPeak( noiseHeights );
  highNoiseHeights = noiseHeights(noiseHeights >= peakNoise) - peakNoise;
  % assume peak is ~ gaussian, and estimate sigma of that peak by finding the
  % location halfway down the cumulative distribution
  numSigmaHalf = sqrt( 2.0 ) * erfinv( 0.5 );
  sigma = median( highNoiseHeights ) / numSigmaHalf;
  numSamplePoints = numel( highNoiseHeights );
  % choose threshold so rare the the probability of a false spike in whole data
  % set is pFalseSpike. do numerically more stable version of this:
  % rareness = 1 - (1 - pFalseSpike).^(1.0 / numel(noiseHeights));
  rareness = -expm1( log1p( -options.pFalseSpike ) ) / numSamplePoints;
  rareness = -expm1( log1p( -options.pFalseSpike ) );
  numSigmaNeeded = sqrt( 2 ) * erfcinv( rareness );
	
  noiseHeight = peakNoise + sigma * numSigmaNeeded;

  if options.debugPlots
    titleStr = makeTitle( 'Spike Thresholds', options );
    
    numPoints = numel( noiseHeights );
    numBins = max(100, round( sqrt( numPoints ) ));
    i1 = 1 + round( (numPoints - 1) * 0.05 ); h1 = noiseHeights(i1);
    i2 = 1 + round( (numPoints - 1) * 0.95 ); h2 = noiseHeights(i2);
    dH = (h2 - h1) / numBins;
    x = 0:dH:max(noiseHeights);
    density = ksdensity( noiseHeights, x );
    %{
    [n, x] = hist(noiseHeights, numBins);
    n = n ./ max(n);
    %}
    fig = NamedFigure(titleStr); fig.WindowStyle = 'docked';
    ax = subplot(1,2,2, 'Parent', fig);
    bar(ax, x, density, 1.0, 'EdgeColor', 'b', 'FaceColor', 'b');
    hold( ax, 'on' )
    plot(ax, [noiseHeight, noiseHeight], [0, 1], 'g')
    hold(ax, 'off')
    xlabel(ax, 'Noise (mV)')
    ylabel(ax, 'Relative Frequency')
    titleStr = makeTitle('Spike Height Threshold', options);
    title(ax, RealUnderscores(titleStr))
    legend(ax, 'Noise', 'Spike height threshold', 'Location', 'Best')
    axis( ax, 'tight' )
    xRange = xlim( ax );
    xRange(2) = min( 3 * noiseHeight, xRange(2) );
    xlim( ax, xRange )
    % we're debugging, so spit out information about the cutoffs
    fprintf( 'GetSpikes.m: spike height cutoff: %g\n', ...
             noiseHeight )
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% high-pass filter a signal
function y = highPassFilter(y, halfFilterLength)
% 1. prepare high-pass filter
filtLen = 1 + 2 * halfFilterLength;
filt = repmat( -1.0 / filtLen, 1, filtLen );
filt(1 + halfFilterLength) = filt(1 + halfFilterLength) + 1.0;
% 2. pad signal symmetrically
y = [flip(y(2:halfFilterLength+1)), ...
     y, ...
     flip(y(end-halfFilterLength-1:end-1))];
% 3. return valid part of convolution between padded-signal and filter
y = conv(y, filt, 'valid');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find the peak of list of data points x
% note - this function assumes x is a sorted row of finite values
function peak = findPeak(x)
  numPts = numel( x );
  i1 = 1 + round( (numPts - 1) * 0.05 );
  i2 = 1 + round( (numPts - 1) * 0.95 );
  numDensityPts = max( 100, round( sqrt( numPts ) ) );
  vals = linspace( x(i1), x(i2), numDensityPts );
  density = ksdensity( x, vals );
  [~, maxInd] = max(density);
  [maxDense, maxInd] = max(density);
  peak = vals(maxInd);  
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
    n1Barrier = n2;
    n1 = n1Barrier;    
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% extend brackets further, to ensure all the features of the spike shape
% can be found
function [n1List, n2List] = extendBrackets( n1List, n2List, v, deriv1, ...
                                            deriv2 )
  n1Barrier = 1; numV = numel( v ); numSpikes = numel( n1List );
  for m = 1:numSpikes
    n1 = n1List(m);
    %while n1 > n1Barrier && ( deriv1(n1) > 0 || v(n1-1) < v(n1) || ...
    %                          deriv2(n1) > 0 )
    while n1 > n1Barrier && deriv1(n1) > 0 && deriv2(n1) > 0
      n1 = n1 - 1;
    end
    while n1 > n1Barrier && deriv2(n1) > max( 0, deriv2(n1-1) )
      n1 = n1 - 1;
    end
    n1List(m) = n1;
    
    n2 = n2List(m);
    if m == numSpikes
      n2Barrier = numV;
    else
      n2Barrier = n1List(m+1);
    end
    %while n2 < n2Barrier && ( deriv1(n2) < 0 || v(n2+1) < v(n2) || ...
    %                          deriv2(n2) > 0 )
    while n2 < n2Barrier && deriv1(n2) > 0 && deriv2(n2) > 0
      n2 = n2 + 1;
    end
    while n2 < n2Barrier && deriv2(n2) > max( 0, deriv2(n2+1) )
      n2 = n2 + 1;
    end
    
    n2List(m) = n2;
    n1Barrier = n2;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function spike = getSpikeShape(n1List, n2List, dT, v, deriv, deriv2, ...
                               options)
minSpikeHeight = getNoiseHeight(v, n1List, n2List, options);
%minSpikeHeight = getNoiseHeightFast(v, n1List, n2List, options);
if options.debugPlots
  fprintf('Noise height = %g\n', minSpikeHeight)
end
minSpikeHeight = max(options.minSpikeHeight, minSpikeHeight);
                             

[n1List, n2List] = extendBrackets( n1List, n2List, v, deriv, deriv2 );
spike.n1List = n1List;
spike.n2List = n2List;

numSpikes = length(n1List);
spike.times = nan(1, numSpikes);
badSpikes = false(1, numSpikes);

%K = deriv2 .* (1 + deriv.^2).^-1.5;
K = deriv2;

spike.maxV.v = nan(1, numSpikes);
spike.maxV.t = nan(1, numSpikes);
spike.maxV.ind = nan(1, numSpikes);
spike.maxDeriv.v = nan(1, numSpikes);
spike.maxDeriv.dV = nan(1, numSpikes);
spike.maxDeriv.t = nan(1, numSpikes);
spike.maxDeriv.ind = nan(1, numSpikes);
spike.minDeriv.v = nan(1, numSpikes);
spike.minDeriv.dV = nan(1, numSpikes);
spike.minDeriv.t = nan(1, numSpikes);
spike.minDeriv.ind = nan(1, numSpikes);
spike.preMinV.v = nan(1, numSpikes);
spike.preMinV.t = nan(1, numSpikes);
spike.preMinV.ind = nan(1, numSpikes);
spike.postMinV.v = nan(1, numSpikes);
spike.postMinV.t = nan(1, numSpikes);
spike.postMinV.ind = nan(1, numSpikes);
spike.preMaxCurve.v = nan(1, numSpikes);
spike.preMaxCurve.K = nan(1, numSpikes);
spike.preMaxCurve.t = nan(1, numSpikes);
spike.preMaxCurve.ind = nan(1, numSpikes);
spike.postMaxCurve.v = nan(1, numSpikes);
spike.postMaxCurve.K = nan(1, numSpikes);
spike.postMaxCurve.t = nan(1, numSpikes);
spike.postMaxCurve.ind = nan(1, numSpikes);
spike.height = nan(1, numSpikes);
spike.width = nan(1, numSpikes);
spike.repolarizationPotential = nan(1, numSpikes);
if numSpikes == 0
  return
end

badSpikeReasons = cell(numSpikes, 1);
for m = 1:numSpikes
  n1 = n1List(m);
  n2 = n2List(m);

  %Find the moment and voltage of maximum depolarization
  [maxV, tMaxV, nMaxV] = getExtremum(v, dT, n1, n2, 'max', false);
  spike.times(m) = tMaxV;
  
  if isnan(tMaxV) || nMaxV == n1 || nMaxV == n2
    badSpikes(m) = true;
    badSpikeReasons{m} = 'Couldn''t bracket spike';
    continue
  end
  
  %Find the max derivative
  [maxDV, tMaxDV, nMaxDV] = ...
    getExtremum(deriv, dT, n1, nMaxV - 1, 'max', true);
  vMaxDV = v(nMaxDV);
  %Find the min derivative
  [minDV, tMinDV, nMinDV] = ...
    getExtremum(deriv, dT, nMaxV + 1, n2, 'min', true);
  vMinDV = v(nMinDV);
  
  %Find the max curvature near the spike
  [preMaxK, tPreMaxK, nPreMaxK] = getExtremum(K, dT, n1, nMaxV-1, ...
					                                    'max', true);
  vPreMaxK = v(nPreMaxK);
  [postMaxK, tPostMaxK, nPostMaxK] = getExtremum(K, dT, nMaxV+1, n2, ...
                                                 'max', true);
  vPostMaxK = v(nPostMaxK);

  %Find minimum voltage before and after spike
  while n1 > 1 && v(n1-1) <= v(n1)
    n1 = n1 - 1;
  end
  while n2 < length(v) && v(n2+1) <= v(n2)
    n2 = n2 + 1;
  end
  [preMinV, tPreMin, nPreMin] = getExtremum(v, dT, n1, n1+3, 'min', true);
  [postMinV, tPostMin, nPostMin] = ...
    getExtremum(v, dT, n2-3, n2, 'min', true);
  
  %height = maxV - min(vPreMaxK, vPostMaxK);
  %height = maxV - vPreMaxK;
  checkHeight = maxV - max(vPreMaxK, vPostMaxK);
  if checkHeight < minSpikeHeight
    % this spike is bad
    badSpikes(m) = true;
    badSpikeReasons{m} = sprintf('spike height too short (%g/%g)', ...
                                 checkHeight, minSpikeHeight);
    continue
  end
  height = maxV - vPreMaxK; % this is the relevant height
  rpp = vPostMaxK - vPreMaxK; % this is the repolarization potential
  
  width = tMinDV - tMaxDV;
  aspect = height / width;
  if aspect < options.minSpikeAspect
    % this spike is bad
    badSpikes(m) = true;
    badSpikeReasons{m} = sprintf('spike is too short and wide (%g/%g)', ...
                                 aspect, options.minSpikeAspect);
  end
  
  spike.maxV.v(m) = maxV;
  spike.maxV.t(m) = tMaxV;
  spike.maxV.ind(m) = nMaxV;
  spike.maxDeriv.v(m) = vMaxDV;
  spike.maxDeriv.dV(m) = maxDV;
  spike.maxDeriv.t(m) = tMaxDV;
  spike.maxDeriv.ind(m) = nMaxDV;
  spike.minDeriv.v(m) = vMinDV;
  spike.minDeriv.dV(m) = minDV;
  spike.minDeriv.t(m) = tMinDV;
  spike.minDeriv.ind(m) = nMinDV;
  spike.preMinV.v(m) = preMinV;
  spike.preMinV.t(m) = tPreMin;
  spike.preMinV.ind(m) = nPreMin;
  spike.postMinV.v(m) = postMinV;
  spike.postMinV.t(m) = tPostMin;
  spike.postMinV.ind(m) = nPostMin;
  spike.preMaxCurve.v(m) = vPreMaxK;
  spike.preMaxCurve.K(m) = preMaxK;
  spike.preMaxCurve.t(m) = tPreMaxK;
  spike.preMaxCurve.ind(m) = nPreMaxK;
  spike.postMaxCurve.v(m) = vPostMaxK;
  spike.postMaxCurve.K(m) = postMaxK;
  spike.postMaxCurve.t(m) = tPostMaxK;
  spike.postMaxCurve.ind(m) = nPostMaxK;
  spike.height(m) = height;
  spike.width(m) = width;
  spike.repolarizationPotential(m) = rpp;
end

if options.removeOutliers
  % first check for extremely short spikes
  spikeHeight = spike.height(~badSpikes);
  medianHeight = median(spikeHeight);
  thresholdHeight = min(0.5 * medianHeight, ...
                        medianHeight - 3 * std(spikeHeight));
  badSpikes = badSpikes | (spike.height < thresholdHeight);
  
  % next check for spikes with very low derivative
  spikeDV = spike.maxDeriv.dV(~badSpikes);
  medianDV = median(spikeDV);
  thresholdDV = min(0.5 * medianDV, medianDV - 3 * std(spikeDV));
  badSpikes = badSpikes | (spike.maxDeriv.dV < thresholdDV);
  if options.debugPlots
    % we're debugging, so print out some information about rejected spikes
    for n = 1:length(badSpikes)
      if badSpikes(n)
        if spike.height(n) < thresholdHeight
          badSpikeReasons{n} = 'short spike height';
        end
        
        badTime = spike.times(n);
        if spike.maxDeriv.dV(n) < thresholdDV
          badSpikeReasons{n} = 'small maxDeriv';
        end
        fprintf('Bad spike at t=%g. Reason %s\n', badTime / 1000, ...
          badSpikeReasons{n})
      end
    end
  end
end

if any(badSpikes)
  % remove bad spikes from spike struct
  
  goodSpikes = ~badSpikes;
  
  fNames1 = fieldnames(spike);
  for n1 = 1:length(fNames1)
    name1 = fNames1{n1};
    try
      fNames2 = fieldnames(spike.(name1));
    catch %#ok<CTCH>
      checkList = spike.(name1);
      spike.(name1) = checkList(goodSpikes);
      continue
    end
    for n2 = 1:length(fNames2)
      name2 = fNames2{n2};
      checkList = spike.(name1).(name2);
      spike.(name1).(name2) = checkList(goodSpikes);
    end
  end
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [maxV, tMax, nMax] = getExtremum(v, dT, n1, n2, extremumStr, ...
                                          simple)
% from a bracketed extremum, find the actual extreme time and value
if nargin < 6
  simple = false;
end

if strcmpi(extremumStr, 'min')
  [maxV, nMax] = min(v(n1:n2));
else
  [maxV, nMax] = max(v(n1:n2));
end
nMax = nMax + n1 - 1;

if simple || nMax == 1 || nMax == length(v)
  tMax = dT * (nMax - 1);
  return
end

%Refine by modeling trace as parabola
n1 = nMax - 1;
n2 = nMax;
n3 = nMax + 1;
t2 = dT * n1;
t3 = dT * n2;
t1 = t2 - dT;

if v(n1) == v(n2)
  if v(n2) == v(n3)
    maxV = v(n2);
    tMax = dT * (n2 - 1);
    return
  else
    tMax = (t1 + t2) / 2;
    coeff = (v(n2) - v(n3)) / ((t2 - tMax)^2 - (t3 - tMax)^2);
  end
elseif v(n2) == v(n3)
  tMax = (t2 + t3) / 2;
  coeff = (v(n2) - v(n1)) / ((t2 - tMax)^2 - (t1 - tMax)^2);
else
  val1 = (v(n2) - v(n1)) / (v(n2) - v(n3));

  b = 2 * (t2 - t1 + val1 * (t3 - t2));
  c = val1 * (t2*t2 - t3*t3) + t1*t1 - t2*t2;

  tMax = -c / b;
  % check for sanity on this extremum time
  if tMax < t1 || t3 < tMax
    tMax = dT * (nMax - 1);
    return
  end

  
  coeff = (v(n2) - v(n1)) / ((t2 - tMax)^2 - (t1 - tMax)^2);
  %arbitrary which formula to use:
  %coeff = (v(n3) - v(n1)) / ((t(n3) - tMax)^2 - (t(n1) - tMax)^2);
end

maxV = v(n2) - coeff * (t2 - tMax)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function freq = getSpikeFrequency(times, tFinal)
if isempty(times) || tFinal == 0
  freq = 0;
  return
end

tHalf = .5 * tFinal;
if isempty(find(times > tHalf, 1))
  %Check if there are no events in the second half of the experiment
  %  if so, presumably it just took a LONG time to settle down, so
  %  label the cell as NOT spiking
  freq = 0;
  return
end

numEvents = length(times);
if numEvents == 1
  freq = 1000 * numEvents / tFinal;
else
  freq = 1000 * (numEvents - 1) / (times(end) - times(1));
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
function fig = plotGetSpikeTimes(dT, v, deriv, deriv2, lowCutoff, highCutoff, ...
                               options)
% Plot the derivatives and thresholds, showing how the affect spike
% detection
titleStr = makeTitle('Derivatives', options);

fig = NamedFigure(titleStr); fig.WindowStyle = 'docked'; clf( fig )
ax = subplot( 2, 1, 1, 'Parent', fig );

numV = length(v);
dTSeconds = 0.001 * dT;
tFinal = dTSeconds * (numV - 1);
plot(ax, 0:dTSeconds:tFinal, deriv, 'b-')
hold(ax, 'on' )
plot(ax, [0, tFinal], [lowCutoff, lowCutoff], 'r-')
plot(ax, [0, tFinal], [highCutoff, highCutoff], 'g-')
%xlabel(ax, 'Time (s)', 'FontSize', 18)
ylabel(ax, 'dV/dT (mV/ms)', 'FontSize', 18)
%title(ax, RealUnderscores(titleStr), 'FontSize', 18)
legend( ax, {'dV/dT', 'low threshold', 'high threshold'}, ...
        'Location', 'NorthOutside', 'Orientation', 'Horizontal' )
hold(ax, 'off') ; axis( ax, 'tight' )
ax.Tag = titleStr;

ax = subplot( 2, 1, 2, 'Parent', fig );
plot( ax, 0:dTSeconds:tFinal, deriv2, 'b-' )
xlabel(ax, 'Time (s)', 'FontSize', 18)
ylabel(ax, 'd^2V/dT^2 (mV/ms^2)', 'FontSize', 18)
axis( ax, 'tight' )
ax.Tag = titleStr;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function titleStr = makeTitle(titleBase, options)
% set the full title for a figure based on base title and plotSubject
if ischar(options.plotSubject)
  titleStr = [options.plotSubject, ': ', titleBase];
else
  titleStr = titleBase;
end
end
