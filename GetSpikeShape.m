%spike = GetSpikeShape( n1List, n2List, dT, v, deriv, deriv2, ...
%                                options )
% From list of bracketed spikes, return spike struct with spike timing and
% shape information
function spike = GetSpikeShape( n1List, n2List, dT, v, deriv, deriv2, ...
                                varargin )
  parser = inputParser();
  parser.KeepUnmatched = true;
  parser.addParameter( 'noiseCheckQuantile', 0.67 )
  parser.addParameter( 'plotSubject', false )
  parser.addParameter( 'debugPlots', false )
  parser.addParameter( 'minSpikeAspect', 0.0 )
  parser.addParameter( 'minSpikeHeight', 0.0 )
  parser.addParameter( 'pFalseSpike', 1e-3 )
  parser.addParameter( 'bracketWidth', [] )
  parser.addParameter( 'removeOutliers', true )
  parser.addParameter( 'outlierFraction', 0.33 )
  parser.addParameter( 'noiseThreshold', [] )
  parser.addParameter( 'checkHeights', [] )
  
  parser.parse( varargin{:} )
  options = parser.Results;
  
  if isempty( deriv )
    % calculate derivatives for spike shape info
    [deriv, deriv2] = getDerivs( n1List, n2List, dT, v, options );
  end
  
  if isempty( options.noiseThreshold ) || isempty( options.checkHeights )
    [minSpikeHeight, checkHeights] = getNoiseHeight( v, n1List, n2List, ...
                                                     options );
  else
    minSpikeHeight = max( options.minSpikeHeight, options.noiseThreshold );
    checkHeights = options.checkHeights;
  end
  %minSpikeHeight = getNoiseHeightFast(v, n1List, n2List, options);
  if options.debugPlots
    fprintf( 'Noise height = %g\n', minSpikeHeight )
  end
  minSpikeHeight = max( options.minSpikeHeight, minSpikeHeight );
  
  
  [n1List, n2List] = extendBrackets( n1List, n2List, v, deriv, deriv2 );
  [spike, numSpikes] = initializeSpike( n1List, n2List );
  if numSpikes == 0
    spike.frequencies = []; spike.intervals = []; spike.freq = 0;
    return
  end
  badSpikes = false(1, numSpikes);
  
  %K = deriv2 .* (1 + deriv.^2).^-1.5;
  K = deriv2;
  badSpikeReasons = cell( numSpikes, 1 );
  for m = 1:numSpikes
    n1 = n1List(m);
    n2 = n2List(m);
    
    %Find the moment and voltage of maximum depolarization
    [maxV, tMaxV, nMaxV] = getExtremum( v, dT, n1, n2, 'max', false );
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
      getExtremum( deriv, dT, nMaxV + 1, n2, 'min', true );
    vMinDV = v(nMinDV);
    
    %Find the max curvature near the spike
    [preMaxK, tPreMaxK, nPreMaxK] = getExtremum( K, dT, n1, nMaxV-1, ...
                                                 'max', true );
    vPreMaxK = v(nPreMaxK);
    [postMaxK, tPostMaxK, nPostMaxK] = getExtremum( K, dT, nMaxV+1, n2, ...
                                                    'max', true);
    vPostMaxK = v(nPostMaxK);
    
    %Find minimum voltage before and after spike
    [preMinV, tPreMin, nPreMin] = ...
      getExtremum( v, dT, n1, n1+3, 'min', true );
    [postMinV, tPostMin, nPostMin] = ...
      getExtremum(v, dT, n2-3, n2, 'min', true);
    
    %height = maxV - min(vPreMaxK, vPostMaxK);
    %height = maxV - vPreMaxK;
    %checkHeight = maxV - max(vPreMaxK, vPostMaxK);
    checkHeight = checkHeights(m);
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
    medianHeight = median( spikeHeight );
    thresholdHeight = options.outlierFraction * medianHeight;
    badSpikes = badSpikes | (spike.height < thresholdHeight);
    
    % next check for spikes with very low derivative
    spikeDV = spike.maxDeriv.dV(~badSpikes);
    medianDV = median( spikeDV );
    thresholdDV = min(0.5 * medianDV, medianDV - 3 * std( spikeDV ));
    badSpikes = badSpikes | (spike.maxDeriv.dV < thresholdDV);
    if options.debugPlots
      % we're debugging, so print out some information about rejected spikes
      for n = 1:numel( badSpikes )
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
  
  if any( badSpikes )
    % remove bad spikes from spike struct
    spike = removeBadSpikes( spike, badSpikes );
  end
  spike.ind = spike.maxV.ind;
  
  %  Calculate spike intervals and frequencies
  if isempty( spike.times )
    spike.intervals = [];
    spike.frequencies = [];
  else
    spike.intervals = spike.times(2:end) - spike.times(1:(end-1));
    spike.frequencies = 1000 ./ spike.intervals;
  end
  %get the overall spike frequency
  spike.freq = getSpikeFrequency( spike.times, dT * (numel( v ) - 1) );
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate derivatives for spike shape info
function [deriv, deriv2] = getDerivs( n1List, n2List, dT, v, options )
  if isempty( options.bracketWidth )
    if isempty( n1List )
      maxTimeWidth = 3.0;
    else
      maxTimeWidth = dT * median( n2List - n1List );
    end
  else
    maxTimeWidth = options.bracketWidth;
  end
  nyquistRate = 1.0 / (2 * dT);
  fStop = min( nyquistRate * 2/3, 1.0 / maxTimeWidth );
  fPass = fStop;
  %nyquistFrac = fStop / nyquistRate;
  [deriv, deriv2] = DerivFilter(v, dT, fPass, fStop);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create structure to hold spike information
function   [spike, numSpikes] = initializeSpike( n1List, n2List )
  spike.n1List = n1List;
  spike.n2List = n2List;
  
  numSpikes = numel( n1List );
  spike.times = nan( 1, numSpikes );
  spike.ind = nan( 1, numSpikes );

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
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% quick method of estimating the height of noise in the trace
function [noiseHeight, checkHeights] = getNoiseHeight(v, n1List, n2List,...
                                                      options)
  if isempty( n1List )
    noiseHeight = 0;
    checkHeights = [];
    return
  end
  
  % find indices when a spike is in progress
  spikeInds = arrayfun( @(n1,n2) n1:n2, n1List, n2List, ...
                        'UniformOutput', false );
  % how wide are spikes
  spikeWidth = median( n2List - n1List );

  % get noise heights as v - lower envelope (v )
  % set filter length to odd integer ~ 1/2 spike width
  filtLen = 1 + 2 * ceil( (spikeWidth - 1) / 8 );
  % get lower-envelope of trace
  envelope = LowerEnvelope( v, filtLen, 'plot', options.debugPlots, ...
                           'title', makeTitle( 'LowerEnvelope', options ) );
  % vFast is the height above envelope
  vFast = v - envelope;
  checkHeights = cellfun( @(inds) max( vFast(inds) ), spikeInds );
  spikeInds = cat(2, spikeInds{:});
  
  minNumDataPoints = 100; %need this many data points to make an okay guess
  fewDataPoints = numel( spikeInds ) > numel( v ) - minNumDataPoints;
  if fewDataPoints
    % not a lot of non-spike data to work with.
    % assume spikes should be more than just large single-point fluctuations,
    % so we want outliers on individual point-to-point differences
    noiseHeights = abs( diff( v ) );
  else % enough spike data
    % spike heights 
    %remove spike indices from consideration
    noiseHeights = vFast;
    noiseHeights(spikeInds) = [];
    
    zeroInds = find( noiseHeights == 0 );
    noiseHeights = arrayfun( @(i1,i2) max( noiseHeights(i1:i2) ), ...
                             [1, zeroInds], ...
                             [zeroInds, numel( noiseHeights )] );
  end

  noiseHeights = noiseHeights(noiseHeights(:) > 0);
  noiseHeights = sort(noiseHeights);
  % find the peak of those noise heights
  [peakNoise, ~, sigma] = FindPeak( noiseHeights, options.noiseCheckQuantile );
  highNoiseHeights = noiseHeights(noiseHeights >= peakNoise) - peakNoise;
  % assume peak is ~ gaussian, and estimate sigma of that peak by finding the
  % location halfway down the cumulative distribution
  
  %numSigmaCheck = sqrt( 2.0 ) * erfinv( options.noiseCheckQuantile );
  %sigma = quantile( highNoiseHeights, options.noiseCheckQuantile ) / numSigmaCheck;
  %numSamplePoints = numel( highNoiseHeights );
  
  % choose threshold so rare the the probability of a false spike in whole data
  % set is pFalseSpike. do numerically more stable version of this:
  % rareness = 1 - (1 - pFalseSpike).^(1.0 / numel(noiseHeights));
  %rareness = -expm1( log1p( -options.pFalseSpike ) ) / numSamplePoints;
  rareness = real( -expm1( log1p( -options.pFalseSpike ) ) );
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

%{
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
      n2Barrier = n1List(m+1) - 1;
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
    n1Barrier = n2 + 1;
  end
end
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [n1List, n2List] = extendBrackets( n1List, n2List, v, foo, bar )
  leftBarrier = 0;
  for spikeInd = 1:numel( n1List )
    n1 = n1List(spikeInd);
    while n1-1 > leftBarrier
      if v(n1-1) < v(n1)
        n1 = n1-1;
      elseif n1-2 > leftBarrier && v(n1-2) < v(n1)
        n1 = n1-2;
      else
        break
      end
    end
    n1List(spikeInd) = n1;
    
    if spikeInd == numel( n1List )
      rightBarrier = numel( v ) + 1;
    else
      rightBarrier = n1List(spikeInd+1);
    end
    n2 = n2List(spikeInd);
    while n2+1 < rightBarrier
      if v(n2+1) < v(n2)
        n2 = n2+1;
      elseif n2+2 < rightBarrier && v(n2+2) < v(n2)
        n2 = n2+2;
      else
        break
      end
    end
    n2List(spikeInd) = n2;
    
    leftBarrier = n2List(spikeInd);
  end
  
  % try to extend n2 past AHP
  for spikeInd = 1:numel( n2List )
    if spikeInd == numel( n1List )
      rightBarrier = numel( v ) + 1;
    else
      rightBarrier = n1List(spikeInd+1);
    end
    n1 = n1List(spikeInd);
    n2 = n2List(spikeInd);
    while n2+1 < rightBarrier
      n2Check = n2 + round( 0.5 * (n2 - n1 ) );
      n2Check = min( max( n2Check, n2+1 ), rightBarrier - 1 );
      [vMin, minInd] = min( v(n2+1:n2Check) );
      if vMin < v(n2)
        n2 = n2 + minInd;
      else
        break
      end
    end
    n2List(spikeInd) = n2;
    
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% from a bracketed extremum, find the actual extreme time and value
function [maxV, tMax, nMax] = getExtremum( v, dT, n1, n2, extremumStr, ...
                                           simple )
  if nargin < 6
    simple = false;
  end

  if strcmpi( extremumStr, 'min' )
    [maxV, nMax] = min( v(n1:n2) );
  else
    [maxV, nMax] = max( v(n1:n2) );
  end
  nMax = nMax + n1 - 1;
  
  if simple || nMax == 1 || nMax == numel( v )
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
% set the full title for a figure based on base title and plotSubject
function titleStr = makeTitle( titleBase, options )
  if ischar( options.plotSubject )
    titleStr = [options.plotSubject, ': ', titleBase];
  else
    titleStr = titleBase;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% remove bad spikes from spike struct
function spike = removeBadSpikes( spike, badSpikes )
  goodSpikes = ~badSpikes;
    
  fNames1 = fieldnames( spike );
  for n1 = 1:numel( fNames1 )
    name1 = fNames1{n1};
    try
      fNames2 = fieldnames( spike.(name1) );
    catch %#ok<CTCH>
      checkList = spike.(name1);
      spike.(name1) = checkList(goodSpikes);
      continue
    end
    for n2 = 1:numel( fNames2 )
      name2 = fNames2{n2};
      checkList = spike.(name1).(name2);
      spike.(name1).(name2) = checkList(goodSpikes);
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function freq = getSpikeFrequency(times, tFinal)
  if isempty( times ) || tFinal == 0
    freq = 0;
    return
  end
  
  tHalf = .5 * tFinal;
  if isempty( find( times > tHalf, 1 ) )
    %Check if there are no events in the second half of the experiment
    %  if so, presumably it just took a LONG time to settle down, so
    %  label the cell as NOT spiking
    freq = 0;
    return
  end
  
  numEvents = numel( times );
  if numEvents == 1
    freq = 1000 * numEvents / tFinal;
  else
    freq = 1000 * (numEvents - 1) / (times(end) - times(1));
  end
end